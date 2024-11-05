# src/utils.jl

"""
    is_axis_aligned(camera::AbstractCamera)

Check if camera pixel grid is aligned with coordinate axes.
Currently assumes IdealCamera is always axis-aligned.

Returns `true` if the camera grid is aligned with the coordinate axes.
"""
function is_axis_aligned(camera::AbstractCamera)
    return true  # For now, assume ideal camera is axis-aligned
    # Later we might want to check for rotated cameras:
    # - Check if x edges vary only in x
    # - Check if y edges vary only in y
end

"""
    get_pixel_size(camera::AbstractCamera)

Get the pixel size from camera geometry.
Assumes uniform pixel size and square pixels.

Returns the pixel size in physical units (typically microns).

# Throws
- AssertionError if pixels are not square
"""
function get_pixel_size(camera::AbstractCamera)
    dx = camera.pixel_edges_x[2] - camera.pixel_edges_x[1]
    dy = camera.pixel_edges_y[2] - camera.pixel_edges_y[1]
    @assert abs(dx - dy) < 1e-10 "Non-square pixels not currently supported"
    return dx
end

"""
    get_pixel_centers(camera::AbstractCamera)

Get arrays of pixel center coordinates from camera geometry.

Returns a tuple of (xcenters, ycenters) where each array contains
the center coordinates of pixels in physical units.
"""
function get_pixel_centers(camera::AbstractCamera)
    xedges = camera.pixel_edges_x
    yedges = camera.pixel_edges_y
    
    xcenters = (xedges[1:end-1] .+ xedges[2:end]) ./ 2
    ycenters = (yedges[1:end-1] .+ yedges[2:end]) ./ 2
    
    return xcenters, ycenters
end

# This one might stay internal as it's mainly for debugging
"""
    _check_normalization(values; tol=1e-6)

Internal function to verify array sums to 1 within tolerance.
"""
function _check_normalization(values; tol=1e-6)
    s = sum(values)
    if abs(s - 1) > tol
        @warn "Values not normalized, sum = $s"
    end
    return s
end

# Export the public functions
export get_pixel_size, get_pixel_centers, is_axis_aligned