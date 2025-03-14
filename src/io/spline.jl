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

"""
    _load_psf_impl(file::HDF5.File, ::Type{SplinePSF})

Load a SplinePSF from an HDF5 file, reconstructing it using the standard constructor.
"""
function _load_psf_impl(file::HDF5.File, ::Type{SplinePSF})
    # Load parameters
    params = file["parameters"]
    _io_check_required_fields(params, ["interp_order", "dimensions"])
    
    interp_order = read(params["interp_order"])
    dimensions = read(params["dimensions"])
    
    # Load coordinate ranges
    x_range = _io_load_range(params, "x_range")
    y_range = _io_load_range(params, "y_range")
    
    # Load the original grid data
    data = file["data"]
    _io_check_required_fields(data, ["original_grid"])
    grid = read(data["original_grid"])
    
    # Create the appropriate SplinePSF based on dimensionality
    if dimensions == "3D"
        # Load z_range for 3D PSF
        z_range = _io_load_range(params, "z_range")
        
        # Create 3D SplinePSF using the standard constructor
        return SplinePSF(grid, x_range, y_range, z_range; order=interp_order)
    else
        # Create 2D SplinePSF
        return SplinePSF(grid, x_range, y_range; order=interp_order)
    end
end