# src/io/io.jl

"""
    IO Module for MicroscopePSFs

Provides unified functionality for saving and loading PSF objects and pupil functions to/from HDF5 files.
All types share a common container format for better interoperability.

Key functions:
- save_psf: Save any PSF or pupil type to an HDF5 file
- load_psf: Load any PSF or pupil type from an HDF5 file
"""

# Import Base functions that we want to extend
import Base: read, write

# Constants
const IO_VERSION = v"1.0.0"  # Version of the I/O system
const MIN_SUPPORTED_VERSION = v"1.0.0"  # Minimum supported version

# Main interface functions
"""
    save_psf(filename::String, psf::AbstractPSF; metadata::Dict=Dict())

Save a PSF to an HDF5 file.

# Arguments
- `filename`: Path where the PSF will be saved
- `psf`: PSF to save
- `metadata`: Optional dictionary of additional metadata to include

# Returns
- `filename` for chaining

# Examples
```julia
# Save a PSF with metadata
psf = Airy2D(1.4, 0.532)
save_psf("airy_psf.h5", psf, metadata=Dict("description" => "Example Airy PSF"))

# Save a PSF with Zernike aberrations
zc = ZernikeCoefficients(15)
add_astigmatism!(zc, 0.5)
psf = Scalar3DPSF(1.4, 0.532, 1.518; coeffs=zc)
save_psf("aberrated_psf.h5", psf)
```
"""
function save_psf(filename::String, psf::AbstractPSF; metadata::Dict=Dict())
    h5open(filename, "w") do file
        # Write common metadata
        attrs = attributes(file)
        attrs["io_version"] = string(IO_VERSION)
        attrs["psf_type"] = string(typeof(psf))
        attrs["creation_date"] = string(Dates.now())
        
        # Write user-provided metadata
        for (key, value) in metadata
            attrs[string(key)] = string(value)
        end
        
        # Dispatch to type-specific save
        _save_psf_impl(file, psf)
    end
    return filename
end

# Generic load interface
"""
    load_psf(filename::String) -> AbstractPSF

Load a PSF from an HDF5 file.

# Arguments
- `filename`: Path to the saved PSF file

# Returns
- PSF of the appropriate type based on stored metadata

# Examples
```julia
# Load a previously saved PSF
psf = load_psf("my_psf.h5")
```
"""
function load_psf(filename::String)
    h5open(filename, "r") do file
        # Check version compatibility
        attrs = attributes(file)
        if haskey(attrs, "io_version")
            _check_version_compatibility(read(attrs["io_version"]))
        end
        
        # Read PSF type
        psf_type_str = read(attrs["psf_type"])
        
        # Convert string to actual type
        try
            psf_type = eval(Meta.parse(psf_type_str))
            println("Loading PSF of type $psf_type_str")
            # Dispatch to appropriate load implementation
            return _load_psf_impl(file, psf_type)
        catch e
            error("Failed to load PSF of type $psf_type_str: $e")
        end
    end
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

# Replace the individual fallback methods with this more general approach

# Generic fallback for loading PSFs with any parameterization
function _load_psf_impl(file::HDF5.File, ::Type{T}) where {T <: AbstractPSF}
    # Get the base type without type parameters
    base_type = T.name.wrapper
    
    # Read the PSF type from file
    attrs = attributes(file)
    psf_type_str = read(attrs["psf_type"])
    
    # Log the type information for debugging
    @debug "Load request for: $T, base type: $base_type, from file: $psf_type_str"
    
    # Try to load using direct type dispatch
    try
        # For each PSF type, we'll use its basic constructor with data from file
        params = file["parameters"]
        
        if base_type <: Gaussian2D
            _io_check_required_fields(params, ["sigma"])
            σ = read(params["sigma"])
            return Gaussian2D(σ)
            
        elseif base_type <: Airy2D
            _io_check_required_fields(params, ["na", "lambda"])
            nₐ = read(params["na"])
            λ = read(params["lambda"])
            return Airy2D(nₐ, λ)
            
        elseif base_type <: Scalar3DPSF
            _io_check_required_fields(params, ["na", "lambda", "n"])
            nₐ = read(params["na"])
            λ = read(params["lambda"])
            n = read(params["n"])
            
            # Load pupil field
            data = file["data"]
            field = _io_load_complex_array(data, "pupil_field")
            
            # Load Zernike coefficients if available
            zernike_coeffs = _io_load_zernike_coeffs(file)
            
            # Create pupil function from field data
            pupil = PupilFunction(nₐ, λ, n, field)
            
            # Create and return the Scalar3DPSF
            return Scalar3DPSF(nₐ, λ, n, pupil, zernike_coeffs)
            
        elseif base_type <: SplinePSF
            _io_check_required_fields(params, ["interp_order", "dimensions"])
            
            # Read dimensions and interpolation order
            dimensions = read(params["dimensions"])
            interp_order = read(params["interp_order"])
            
            # Load coordinate ranges
            x_range = _io_load_range(params, "x_range")
            y_range = _io_load_range(params, "y_range")
            
            # Load the original grid data
            data = file["data"]
            grid = read(data["original_grid"])
            
            # Create appropriate SplinePSF based on dimensionality
            if dimensions == "3D"
                # Load z_range for 3D PSF
                z_range = _io_load_range(params, "z_range")
                
                # Create 3D SplinePSF
                return SplinePSF(grid, x_range, y_range, z_range; order=interp_order)
            else
                # Create 2D SplinePSF
                return SplinePSF(grid, x_range, y_range; order=interp_order)
            end
            
        elseif base_type <: ZernikeCoefficients
            # Load magnitude and phase arrays
            data = file["data"]
            _io_check_required_fields(data, ["magnitude", "phase"])
            
            mag = read(data["magnitude"])
            phase = read(data["phase"])
            
            # Create and return the ZernikeCoefficients
            return ZernikeCoefficients(mag, phase)
            
        else
            # If no specific implementation exists, use the default error
            error("No implementation for loading $base_type")
        end
    catch e
        @error "Error loading PSF: $e"
        rethrow(e)
    end
end

# Default error implementation
function _save_psf_impl(file::HDF5.File, psf::AbstractPSF)
    error("Saving not implemented for PSF type $(typeof(psf))")
end

function _load_psf_impl(file::HDF5.File, ::Type{T}) where {T <: AbstractPSF}
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