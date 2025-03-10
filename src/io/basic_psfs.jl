# src/io/basic_psfs.jl

"""
I/O implementations for basic PSF types (Gaussian2D and Airy2D).
These are simpler PSFs with minimal parameters.
"""

# ---- Gaussian2D ----

"""
    _save_psf_impl(file::HDF5.File, psf::Gaussian2D)

Save a Gaussian2D PSF to an HDF5 file.
"""
function _save_psf_impl(file::HDF5.File, psf::Gaussian2D)
    params = create_group(file, "parameters")
    params["sigma"] = psf.σ
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{Gaussian2D})

Load a Gaussian2D PSF from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{Gaussian2D})
    params = file["parameters"]
    _io_check_required_fields(params, ["sigma"])
    
    σ = read(params["sigma"])
    return Gaussian2D(σ)
end

# ---- Airy2D ----

"""
    _save_psf_impl(file::HDF5.File, psf::Airy2D)

Save an Airy2D PSF to an HDF5 file.
"""
function _save_psf_impl(file::HDF5.File, psf::Airy2D)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    # Note: No need to save ν as it's derived from nₐ and λ
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{Airy2D})

Load an Airy2D PSF from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{Airy2D})
    params = file["parameters"]
    _io_check_required_fields(params, ["na", "lambda"])
    
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    return Airy2D(nₐ, λ)
end