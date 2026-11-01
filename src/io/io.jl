# src/io/io.jl

"""
    IO Module for MicroscopePSFs

Provides unified functionality for saving and loading PSF objects and pupil functions to/from HDF5 files.
All types share a common container format for better interoperability.

Key functions:
- save_psf: Save any PSF or pupil type to an HDF5 file
- load_psf: Load any PSF or pupil type from an HDF5 file
"""

# Constants
const _H5Parent = Union{HDF5.File, HDF5.Group}
const IO_VERSION = v"1.1.0"  # Version of the I/O system
const MIN_SUPPORTED_VERSION = v"1.0.0"  # Minimum supported version

# Main interface functions
"""
    save_psf(filename::String, object; metadata::Dict=Dict())

Save a PSF or related object (e.g., ZernikeCoefficients, PupilFunction) to an HDF5 file.

# Arguments
- `filename`: Path where the PSF will be saved
- `object`: Object to save (PSF, ZernikeCoefficients, PupilFunction, etc.)
- `metadata`: Optional dictionary of additional metadata to include

# Returns
- `filename` for chaining
"""
function save_psf(filename::String, object; metadata::Dict=Dict())
    h5open(filename, "w") do file
        # Write common metadata
        attrs = attributes(file)
        attrs["io_version"] = string(IO_VERSION)
        
        # Store the bare type name without module or parameters. `string(typeof(object))`
        # would qualify it as "MicroscopePSFs.X" whenever Main has not run
        # `using MicroscopePSFs`.
        attrs["psf_type"] = String(nameof(typeof(object)))
        
        attrs["creation_date"] = string(Dates.now())
        
        # Write user-provided metadata
        for (key, value) in metadata
            attrs[string(key)] = string(value)
        end
        
        # Dispatch to type-specific save
        _save_psf_impl(file, object)
    end
    return filename
end

# Generic load interface
"""
    load_psf(filename::String)

Load a PSF or related object (e.g., ZernikeCoefficients, PupilFunction) from an HDF5 file.

# Arguments
- `filename`: Path to the HDF5 file containing a saved PSF object

# Returns
- The loaded object with its original type (PSF, ZernikeCoefficients, PupilFunction, etc.)

# Examples
```julia
# Load a previously saved PSF
psf = load_psf("my_airy_psf.h5")

# Use the loaded PSF normally
intensity = psf(0.1, 0.2)

# Load Zernike coefficients and use them
zc = load_psf("zernike_coeffs.h5")
psf = ScalarPSF(1.4, 0.532, 1.518; zernike_coeffs=zc)
```

See also: [`save_psf`](@ref)
"""
function load_psf(filename::String)
    h5open(filename, "r") do file
        # Check version compatibility
        attrs = attributes(file)
        file_version = haskey(attrs, "io_version") ?
            read(attrs["io_version"]) : string(MIN_SUPPORTED_VERSION)
        _check_version_compatibility(file_version)
        
        # Read object type as a base type name
        type_str = read(attrs["psf_type"])
        @info "Loading PSF of type $(type_str)"
        
        object_type = _resolve_psf_type(type_str)

        try
            # Dispatch to appropriate load implementation
            result = _load_psf_impl(file, object_type)
            
            return result
        catch e
            error("Failed to load object of type $type_str: $e")
        end
    end
end

"""
    _resolve_psf_type(type_str::AbstractString) -> Type

Map a stored `psf_type` attribute to a MicroscopePSFs type. Accepts the bare name
(`"SplinePSF"`) and the `MicroscopePSFs.`-qualified name found in some saved files.
Any other name throws an `ArgumentError`.
"""
function _resolve_psf_type(type_str::AbstractString)
    name = chopprefix(type_str, "MicroscopePSFs.")
    sym = Symbol(name)
    if occursin('.', name) || !isdefined(MicroscopePSFs, sym) ||
        !(getfield(MicroscopePSFs, sym) isa Type)
        throw(ArgumentError(
            "Unrecognized psf_type \"$type_str\" in PSF file: " *
            "expected a MicroscopePSFs type name such as \"SplinePSF\""))
    end
    return getfield(MicroscopePSFs, sym)
end

# Version compatibility check
function _check_version_compatibility(file_version::String)
    file_ver = VersionNumber(file_version)
    
    if file_ver < MIN_SUPPORTED_VERSION
        error("PSF file version ($file_ver) is too old. Minimum supported version is $(MIN_SUPPORTED_VERSION).")
    elseif file_ver.major > IO_VERSION.major
        @warn "PSF file version ($file_ver) is newer than supported version ($(IO_VERSION)). Some features may not be supported."
    end
end

# Default error implementation
function _save_psf_impl(file::_H5Parent, psf::AbstractPSF)
    error("Saving not implemented for PSF type $(typeof(psf))")
end

function _load_psf_impl(file::_H5Parent, ::Type{T}) where {T <: AbstractPSF}
    error("Loading not implemented for PSF type $T")
end



# Helper files with shared functionality will be included here
include("common.jl")

# Type-specific implementations
include("basic_psfs.jl")
include("scalar3d.jl")
include("spline.jl")
include("vector3d.jl")
include("pupil.jl")
include("zernike.jl")
include("multichannel.jl")
