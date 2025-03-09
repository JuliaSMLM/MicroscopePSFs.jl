# src/io/common.jl

"""
Shared helper functions for PSF I/O operations.
These simplify common patterns when saving and loading different PSF types.
"""

# ------ Complex array helpers ------

"""
    _io_save_complex_array(group, name, array)

Save a complex array to an HDF5 group by splitting into real and imaginary parts.
"""
function _io_save_complex_array(group, name, array)
    group["$(name)_real"] = real(array)
    group["$(name)_imag"] = imag(array)
end

"""
    _io_load_complex_array(group, name)

Load a complex array from an HDF5 group by combining real and imaginary parts.
"""
function _io_load_complex_array(group, name)
    real_part = read(group["$(name)_real"])
    imag_part = read(group["$(name)_imag"])
    return real_part + im * imag_part
end

# ------ Range helpers ------

"""
    _io_save_range(group, name, range)

Save a range to an HDF5 group by storing start, step, and length.
"""
function _io_save_range(group, name, range)
    group["$(name)_start"] = first(range)
    group["$(name)_step"] = step(range)
    group["$(name)_length"] = length(range)
end

"""
    _io_load_range(group, name)

Load a range from an HDF5 group by reconstructing from start, step, and length.
"""
function _io_load_range(group, name)
    start = read(group["$(name)_start"])
    step_val = read(group["$(name)_step"]) 
    len = read(group["$(name)_length"])
    return range(start, step=step_val, length=len)
end

# ------ Zernike coefficient helpers ------

"""
    _io_save_zernike_coeffs(file, coeffs; group_name="zernike_coefficients")

Save ZernikeCoefficients to an HDF5 file in a specified group.
Does nothing if coeffs is nothing.
"""
function _io_save_zernike_coeffs(file, coeffs; group_name="zernike_coefficients")
    isnothing(coeffs) && return
    
    zernike = create_group(file, group_name)
    zernike["magnitude"] = coeffs.mag
    zernike["phase"] = coeffs.phase
end

"""
    _io_load_zernike_coeffs(file; group_name="zernike_coefficients")

Load ZernikeCoefficients from an HDF5 file from a specified group.
Returns nothing if the group doesn't exist.
"""
function _io_load_zernike_coeffs(file; group_name="zernike_coefficients")
    !haskey(file, group_name) && return nothing
    
    zc = file[group_name]
    mag = read(zc["magnitude"])
    phase = read(zc["phase"])
    return ZernikeCoefficients(mag, phase)
end

# ------ Pupil parameter helpers ------

"""
    _io_save_pupil_params(params, pupil)

Save common pupil parameters to an HDF5 group.
"""
function _io_save_pupil_params(params, pupil)
    params["na"] = pupil.nₐ
    params["lambda"] = pupil.λ
    params["n"] = pupil.n
end

"""
    _io_check_required_fields(group, fields)

Check if all required fields exist in an HDF5 group.
Throws an error if any required field is missing.
"""
function _io_check_required_fields(group, fields)
    for field in fields
        if !haskey(group, field)
            error("Required field '$field' missing from HDF5 file")
        end
    end
end