# src/io.jl

"""
I/O functionality for PSFs including serialization to/from HDF5 files
"""

# Import Base functions that we want to extend
import Base: read, write
import Base: save, load

# Constants
const IO_VERSION = v"1.0.0"  # Version of the I/O system

# Generic save interface
"""
    save_psf(filename::String, psf::AbstractPSF; metadata::Dict=Dict())

Save a PSF to an HDF5 file.

# Arguments
- `filename`: Path where the PSF will be saved
- `psf`: PSF to save
- `metadata`: Optional dictionary of additional metadata to include

# Returns
- `filename` for chaining
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
- PSF of the appropriate type
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
    
    if file_ver.major > IO_VERSION.major
        @warn "PSF file version ($file_ver) is newer than supported version ($(IO_VERSION)). Some features may not be supported."
    end
end

# Default error implementation
function _save_psf_impl(file::HDF5.File, psf::AbstractPSF)
    error("Saving not implemented for PSF type $(typeof(psf))")
end

function _load_psf_impl(file::HDF5.File, ::Type{T}) where {T <: AbstractPSF}
    error("Loading not implemented for PSF type $T")
end

# Implementations for specific PSF types
# Gaussian2D
function _save_psf_impl(file::HDF5.File, psf::Gaussian2D)
    params = create_group(file, "parameters")
    params["sigma"] = psf.σ
end

function _load_psf_impl(file::HDF5.File, ::Type{Gaussian2D})
    params = file["parameters"]
    σ = read(params["sigma"])
    return Gaussian2D(σ)
end

# Airy2D
function _save_psf_impl(file::HDF5.File, psf::Airy2D)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
end

function _load_psf_impl(file::HDF5.File, ::Type{Airy2D})
    params = file["parameters"]
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    return Airy2D(nₐ, λ)
end

# Scalar3DPSF
function _save_psf_impl(file::HDF5.File, psf::Scalar3DPSF)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    params["n"] = psf.n
    
    # Save pupil field
    data = create_group(file, "data")
    data["pupil_field_real"] = real(psf.pupil.field)
    data["pupil_field_imag"] = imag(psf.pupil.field)
end

function _load_psf_impl(file::HDF5.File, ::Type{Scalar3DPSF})
    params = file["parameters"]
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n = read(params["n"])
    
    data = file["data"]
    real_part = read(data["pupil_field_real"])
    imag_part = read(data["pupil_field_imag"])
    field = real_part + im * imag_part
    
    pupil = PupilFunction(nₐ, λ, n, field)
    return Scalar3DPSF(nₐ, λ, n, pupil)
end

# SplinePSF
function _save_psf_impl(file::HDF5.File, psf::SplinePSF)
    params = create_group(file, "parameters")
    
    # Save range information
    params["x_start"] = first(psf.x_range)
    params["x_step"] = step(psf.x_range)
    params["x_length"] = length(psf.x_range)
    
    params["y_start"] = first(psf.y_range)
    params["y_step"] = step(psf.y_range)
    params["y_length"] = length(psf.y_range)
    
    # Handle 2D vs 3D
    if psf.z_range !== nothing
        params["dimensions"] = "3D"
        params["z_start"] = first(psf.z_range)
        params["z_step"] = step(psf.z_range)
        params["z_length"] = length(psf.z_range)
        
        # Sample the PSF on the grid
        grid = Array{Float64}(undef, length(psf.y_range), length(psf.x_range), length(psf.z_range))
        for (iz, z) in enumerate(psf.z_range)
            for (ix, x) in enumerate(psf.x_range)
                for (iy, y) in enumerate(psf.y_range)
                    grid[iy, ix, iz] = psf(x, y, z)
                end
            end
        end
    else
        params["dimensions"] = "2D"
        
        # Sample the PSF on the grid
        grid = Array{Float64}(undef, length(psf.y_range), length(psf.x_range))
        for (ix, x) in enumerate(psf.x_range)
            for (iy, y) in enumerate(psf.y_range)
                grid[iy, ix] = psf(x, y)
            end
        end
    end
    
    # Save the sampled grid
    data = create_group(file, "data")
    data["spline_grid"] = grid
end

function _load_psf_impl(file::HDF5.File, ::Type{SplinePSF})
    params = file["parameters"]
    
    # Read range information
    x_start = read(params["x_start"])
    x_step = read(params["x_step"])
    x_length = read(params["x_length"])
    
    y_start = read(params["y_start"])
    y_step = read(params["y_step"])
    y_length = read(params["y_length"])
    
    # Create ranges
    x_range = range(x_start, step=x_step, length=x_length)
    y_range = range(y_start, step=y_step, length=y_length)
    
    # Read grid data
    data = file["data"]
    grid = read(data["spline_grid"])
    
    # Handle 2D vs 3D
    dimensions = read(params["dimensions"])
    if dimensions == "3D"
        z_start = read(params["z_start"])
        z_step = read(params["z_step"])
        z_length = read(params["z_length"])
        z_range = range(z_start, step=z_step, length=z_length)
        
        return SplinePSF(grid, x_range, y_range, z_range)
    else
        return SplinePSF(grid, x_range, y_range)
    end
end

# Vector3DPSF
function _save_psf_impl(file::HDF5.File, psf::Vector3DPSF)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    params["n_medium"] = psf.n_medium
    params["n_coverslip"] = psf.n_coverslip
    params["n_immersion"] = psf.n_immersion
    params["focal_z"] = psf.focal_z
    
    # Save dipole orientation
    dipole = create_group(file, "dipole")
    dipole["px"] = psf.dipole.px
    dipole["py"] = psf.dipole.py
    dipole["pz"] = psf.dipole.pz
    
    # Save pupil functions
    data = create_group(file, "data")
    pupil_ex = create_group(data, "pupil_ex")
    pupil_ex["field_real"] = real(psf.pupil.Ex.field)
    pupil_ex["field_imag"] = imag(psf.pupil.Ex.field)
    
    pupil_ey = create_group(data, "pupil_ey")
    pupil_ey["field_real"] = real(psf.pupil.Ey.field)
    pupil_ey["field_imag"] = imag(psf.pupil.Ey.field)
end

function _load_psf_impl(file::HDF5.File, ::Type{Vector3DPSF})
    params = file["parameters"]
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n_medium = read(params["n_medium"])
    n_coverslip = read(params["n_coverslip"])
    n_immersion = read(params["n_immersion"])
    focal_z = read(params["focal_z"])
    
    # Load dipole orientation
    dipole_group = file["dipole"]
    px = read(dipole_group["px"])
    py = read(dipole_group["py"])
    pz = read(dipole_group["pz"])
    dipole = DipoleVector(px, py, pz)
    
    # Load pupil fields
    data = file["data"]
    pupil_ex = data["pupil_ex"]
    ex_real = read(pupil_ex["field_real"])
    ex_imag = read(pupil_ex["field_imag"])
    ex_field = ex_real + im * ex_imag
    
    pupil_ey = data["pupil_ey"]
    ey_real = read(pupil_ey["field_real"])
    ey_imag = read(pupil_ey["field_imag"])
    ey_field = ey_real + im * ey_imag
    
    # Create pupil functions
    Ex = PupilFunction(nₐ, λ, n_medium, ex_field)
    Ey = PupilFunction(nₐ, λ, n_medium, ey_field)
    vpupil = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, Ex, Ey)
    
    return Vector3DPSF(nₐ, λ, n_medium, n_coverslip, n_immersion, dipole, focal_z, vpupil)
end

# Convenience aliases
"""
    save(filename::String, psf::AbstractPSF; kwargs...)

Alias for save_psf.
"""
function save(filename::String, psf::AbstractPSF; kwargs...)
    return save_psf(filename, psf; kwargs...)
end

"""
    load(filename::String)

Alias for load_psf.
"""
function load(filename::String)
    return load_psf(filename)
end

# PupilFunction I/O
"""
    _save_psf_impl(file::HDF5.File, pupil::PupilFunction)

Save a PupilFunction to an HDF5 file.
"""
function _save_psf_impl(file::HDF5.File, pupil::PupilFunction)
    params = create_group(file, "parameters")
    params["na"] = pupil.nₐ
    params["lambda"] = pupil.λ
    params["n"] = pupil.n
    
    # Save field data
    data = create_group(file, "data")
    data["field_real"] = real(pupil.field)
    data["field_imag"] = imag(pupil.field)
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{PupilFunction})

Load a PupilFunction from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{PupilFunction})
    params = file["parameters"]
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n = read(params["n"])
    
    data = file["data"]
    real_part = read(data["field_real"])
    imag_part = read(data["field_imag"])
    field = real_part + im * imag_part
    
    return PupilFunction(nₐ, λ, n, field)
end

# Export the main public functions
export save_psf, load_psf, save, load