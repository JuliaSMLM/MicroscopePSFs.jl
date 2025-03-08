# src/io.jl

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

# ======= Implementations for specific PSF types =======

# ---- Gaussian2D ----
function _save_psf_impl(file::HDF5.File, psf::Gaussian2D)
    params = create_group(file, "parameters")
    params["sigma"] = psf.σ
end

function _load_psf_impl(file::HDF5.File, ::Type{Gaussian2D})
    params = file["parameters"]
    σ = read(params["sigma"])
    return Gaussian2D(σ)
end

# ---- Airy2D ----
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

# ---- Scalar3DPSF ----
function _save_psf_impl(file::HDF5.File, psf::Scalar3DPSF)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    params["n"] = psf.n
    
    # Save Zernike coefficients if available
    if !isnothing(psf.zernike_coeffs)
        zernike = create_group(file, "zernike_coefficients")
        zernike["magnitude"] = psf.zernike_coeffs.mag
        zernike["phase"] = psf.zernike_coeffs.phase
    end
    
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
    
    # Load pupil field
    data = file["data"]
    real_part = read(data["pupil_field_real"])
    imag_part = read(data["pupil_field_imag"])
    field = real_part + im * imag_part
    
    # Load Zernike coefficients if available
    zernike_coeffs = nothing
    if haskey(file, "zernike_coefficients")
        zc = file["zernike_coefficients"]
        mag = read(zc["magnitude"])
        phase = read(zc["phase"])
        zernike_coeffs = ZernikeCoefficients(mag, phase)
    end
    
    pupil = PupilFunction(nₐ, λ, n, field)
    return Scalar3DPSF(nₐ, λ, n, pupil, zernike_coeffs)
end

# ---- SplinePSF ----
function _save_psf_impl(file::HDF5.File, psf::SplinePSF)
    params = create_group(file, "parameters")
    
    # Save range information
    params["x_start"] = first(psf.x_range)
    params["x_step"] = step(psf.x_range)
    params["x_length"] = length(psf.x_range)
    
    params["y_start"] = first(psf.y_range)
    params["y_step"] = step(psf.y_range)
    params["y_length"] = length(psf.y_range)
    
    # Save interpolation order directly
    params["interp_order"] = psf.interp_order
    
    # Handle 2D vs 3D
    if psf.z_range !== nothing
        params["dimensions"] = "3D"
        params["z_start"] = first(psf.z_range)
        params["z_step"] = step(psf.z_range)
        params["z_length"] = length(psf.z_range)
    else
        params["dimensions"] = "2D"
    end
    
    # Save the original grid - this is the key improvement
    data = create_group(file, "data")
    data["original_grid"] = psf.original_grid
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
    
    # Read interpolation order
    interp_order = 3  # Default to cubic
    if haskey(params, "interp_order")
        interp_order = read(params["interp_order"])
    end
    
    # Read the original grid
    data = file["data"]
    grid = read(data["original_grid"])
    
    # Handle 2D vs 3D
    dimensions = read(params["dimensions"])
    if dimensions == "3D"
        z_start = read(params["z_start"])
        z_step = read(params["z_step"])
        z_length = read(params["z_length"])
        z_range = range(z_start, step=z_step, length=z_length)
        
        # Create the SplinePSF with the original grid and interpolation order
        return SplinePSF(grid, x_range, y_range, z_range; order=interp_order)
    else
        # Create the SplinePSF with the original grid and interpolation order
        return SplinePSF(grid, x_range, y_range; order=interp_order)
    end
end

# ---- Vector3DPSF ----
function _save_psf_impl(file::HDF5.File, psf::Vector3DPSF)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    params["n_medium"] = psf.n_medium
    params["n_coverslip"] = psf.n_coverslip
    params["n_immersion"] = psf.n_immersion
    params["focal_z"] = psf.focal_z
    
    # Save Zernike coefficients if available
    if !isnothing(psf.zernike_coeffs)
        zernike = create_group(file, "zernike_coefficients")
        zernike["magnitude"] = psf.zernike_coeffs.mag
        zernike["phase"] = psf.zernike_coeffs.phase
    end
    
    # Save dipole orientation
    dipole = create_group(file, "dipole")
    dipole["px"] = psf.dipole.px
    dipole["py"] = psf.dipole.py
    dipole["pz"] = psf.dipole.pz
    
    # Save pupil functions
    data = create_group(file, "data")
    pupil_ex = create_group(data, "pupil_ex")
    pupil_ex["field_real"] = real(psf.vector_pupils.Ex.field)
    pupil_ex["field_imag"] = imag(psf.vector_pupils.Ex.field)
    
    pupil_ey = create_group(data, "pupil_ey")
    pupil_ey["field_real"] = real(psf.vector_pupils.Ey.field)
    pupil_ey["field_imag"] = imag(psf.vector_pupils.Ey.field)
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
    
    # Load Zernike coefficients if available
    zernike_coeffs = nothing
    if haskey(file, "zernike_coefficients")
        zc = file["zernike_coefficients"]
        mag = read(zc["magnitude"])
        phase = read(zc["phase"])
        zernike_coeffs = ZernikeCoefficients(mag, phase)
    end
    
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
    
    # Create the Vector3DPSF with proper parameters
    # For the base_pupil, we'd need to decide on a reconstruction strategy
    base_pupil = nothing
    
    return Vector3DPSF{typeof(nₐ)}(
        nₐ, λ, n_medium, n_coverslip, n_immersion,
        dipole, focal_z, vpupil, base_pupil, zernike_coeffs
    )
end

# ---- PupilFunction ----
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

# ---- VectorPupilFunction ----
function _save_psf_impl(file::HDF5.File, pupil::VectorPupilFunction)
    params = create_group(file, "parameters")
    params["na"] = pupil.nₐ
    params["lambda"] = pupil.λ
    params["n_medium"] = pupil.n_medium
    params["n_coverslip"] = pupil.n_coverslip
    params["n_immersion"] = pupil.n_immersion
    
    # Save Ex field data
    ex_data = create_group(file, "ex_field")
    ex_data["field_real"] = real(pupil.Ex.field)
    ex_data["field_imag"] = imag(pupil.Ex.field)
    
    # Save Ey field data
    ey_data = create_group(file, "ey_field")
    ey_data["field_real"] = real(pupil.Ey.field)
    ey_data["field_imag"] = imag(pupil.Ey.field)
end

function _load_psf_impl(file::HDF5.File, ::Type{VectorPupilFunction})
    params = file["parameters"]
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n_medium = read(params["n_medium"])
    n_coverslip = read(params["n_coverslip"])
    n_immersion = read(params["n_immersion"])
    
    # Load Ex field
    ex_data = file["ex_field"]
    ex_real = read(ex_data["field_real"])
    ex_imag = read(ex_data["field_imag"])
    ex_field = ex_real + im * ex_imag
    
    # Load Ey field
    ey_data = file["ey_field"]
    ey_real = read(ey_data["field_real"])
    ey_imag = read(ey_data["field_imag"])
    ey_field = ey_real + im * ey_imag
    
    # Create pupil functions
    Ex = PupilFunction(nₐ, λ, n_medium, ex_field)
    Ey = PupilFunction(nₐ, λ, n_medium, ey_field)
    
    return VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, Ex, Ey)
end

# ---- ZernikeCoefficients ----
function _save_psf_impl(file::HDF5.File, zc::ZernikeCoefficients)
    data = create_group(file, "data")
    data["magnitude"] = zc.mag
    data["phase"] = zc.phase
end

function _load_psf_impl(file::HDF5.File, ::Type{ZernikeCoefficients})
    data = file["data"]
    mag = read(data["magnitude"])
    phase = read(data["phase"])
    return ZernikeCoefficients(mag, phase)
end

