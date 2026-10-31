# src/io/basic_psfs.jl

"""
I/O implementations for basic PSF types (GaussianPSF and AiryPSF).
These are simpler PSFs with minimal parameters.
"""

# ---- GaussianPSF ----

"""
    _save_psf_impl(file::_H5Parent, psf::GaussianPSF)

Save a GaussianPSF to an HDF5 file.
"""
function _save_psf_impl(file::_H5Parent, psf::GaussianPSF)
    params = create_group(file, "parameters")
    params["sigma"] = psf.σ
end

"""
    _load_psf_impl(file::_H5Parent, ::Type{GaussianPSF})

Load a GaussianPSF from an HDF5 file.
"""
function _load_psf_impl(file::_H5Parent, ::Type{GaussianPSF})
    params = file["parameters"]
    _io_check_required_fields(params, ["sigma"])
    
    σ = read(params["sigma"])
    return GaussianPSF(σ)
end

# ---- AiryPSF ----

"""
    _save_psf_impl(file::_H5Parent, psf::AiryPSF)

Save an AiryPSF to an HDF5 file.
"""
function _save_psf_impl(file::_H5Parent, psf::AiryPSF)
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    # Note: No need to save ν as it's derived from nₐ and λ
end

"""
    _load_psf_impl(file::_H5Parent, ::Type{AiryPSF})

Load an AiryPSF from an HDF5 file.
"""
function _load_psf_impl(file::_H5Parent, ::Type{AiryPSF})
    params = file["parameters"]
    _io_check_required_fields(params, ["na", "lambda"])
    
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    return AiryPSF(nₐ, λ)
end