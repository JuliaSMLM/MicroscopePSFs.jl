# src/io/spline.jl

"""
I/O implementations for SplinePSF type.
Ensures proper saving and loading of the original data grid
and interpolation parameters needed to reconstruct the spline interpolant.
"""

"""
    _save_psf_impl(file::HDF5.File, psf::SplinePSF)

Save a SplinePSF to an HDF5 file, including the original grid data
used to create the spline interpolation.
"""
function _save_psf_impl(file::HDF5.File, psf::SplinePSF)
    # Store parameters
    params = create_group(file, "parameters")
    
    # Save interpolation order
    params["interp_order"] = psf.interp_order
    
    # Determine dimensionality
    is_3d = psf.z_range !== nothing
    params["dimensions"] = is_3d ? "3D" : "2D"
    
    # Save range information using helpers
    _io_save_range(params, "x_range", psf.x_range)
    _io_save_range(params, "y_range", psf.y_range)
    
    # Save z_range if this is a 3D PSF
    if is_3d
        _io_save_range(params, "z_range", psf.z_range)
    end
    
    # Save the original grid data
    data = create_group(file, "data")
    data["original_grid"] = psf.original_grid
end

# Note: The loading of SplinePSF is handled by the generic _load_psf_impl in io.jl
# We don't need a specific loading function here since:
# 1. The SplinePSF constructor takes the original grid and ranges directly
# 2. The generic implementation in io.jl handles the type dispatch based on the base type
# 3. This approach avoids any dependency on implementation details like OffsetArrays