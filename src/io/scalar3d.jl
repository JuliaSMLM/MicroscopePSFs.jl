# src/io/scalar3d.jl

"""
I/O implementations for ScalarPSF type.
Handles saving and loading of complex pupil field data and Zernike coefficients.
"""

"""
    _save_psf_impl(file::HDF5.File, psf::ScalarPSF)

Save a ScalarPSF to an HDF5 file, including pupil function and any Zernike coefficients.
"""
function _save_psf_impl(file::HDF5.File, psf::ScalarPSF)
    # Save basic parameters
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    params["n"] = psf.n
    
    # Save Zernike coefficients if available
    if !isnothing(psf.zernike_coeffs)
        _io_save_zernike_coeffs(file, psf.zernike_coeffs)
    end
    
    # Save pupil field
    data = create_group(file, "data")
    _io_save_complex_array(data, "pupil_field", psf.pupil.field)
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{ScalarPSF}) 

Load a ScalarPSF from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{ScalarPSF}) 
    # Load parameters
    params = file["parameters"]
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
    
    # Create and return the ScalarPSF
    return ScalarPSF(nₐ, λ, n, pupil, zernike_coeffs)
end