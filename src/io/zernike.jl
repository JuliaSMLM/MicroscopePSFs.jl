# src/io/zernike.jl

"""
I/O implementations for ZernikeCoefficients type.
Enables saving and loading of Zernike polynomial coefficients
used to represent pupil aberrations.
"""

"""
    _save_psf_impl(file::HDF5.File, zc::ZernikeCoefficients)

Save ZernikeCoefficients to an HDF5 file, storing magnitude and phase arrays.
"""
function _save_psf_impl(file::HDF5.File, zc::ZernikeCoefficients)
    # Create data group
    data = create_group(file, "data")
    
    # Save magnitude and phase arrays
    data["magnitude"] = zc.mag
    data["phase"] = zc.phase
    
    # Save indexing information
    params = create_group(file, "parameters")
    params["num_terms"] = length(zc.mag)
    params["max_radial_order"] = max_radial_order(length(zc.mag))
    
    # Store non-zero terms for quick reference
    significant = significant_terms(zc)
    if !isempty(significant)
        sig_group = create_group(data, "significant")
        sig_group["indices"] = [t[1] for t in significant]
        sig_group["mag_values"] = [t[2] for t in significant]
        sig_group["phase_values"] = [t[3] for t in significant]
    end
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{ZernikeCoefficients{T}}) where T

Load ZernikeCoefficients from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{ZernikeCoefficients{T}}) where T
    # Load magnitude and phase arrays
    data = file["data"]
    _io_check_required_fields(data, ["magnitude", "phase"])
    
    mag = read(data["magnitude"])
    phase = read(data["phase"])
    
    # Create and return the ZernikeCoefficients
    return ZernikeCoefficients(mag, phase)
end