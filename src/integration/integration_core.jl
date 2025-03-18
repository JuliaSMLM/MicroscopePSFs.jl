# src/integration/integration_core.jl

#=
PSF and Emitter Trait System
This allows automatic detection of PSF capabilities and emitter properties
=#

"""
    supports_3d(psf_type::Type{<:AbstractPSF}) -> Bool

Check if a PSF type supports 3D evaluation via psf(x, y, z).
"""
supports_3d(::Type{<:AbstractPSF}) = false
supports_3d(::Type{<:Abstract3DPSF}) = true

# Special case for SplinePSF - depends on whether z_range is defined
supports_3d(::Type{SplinePSF{T,IT}}) where {T,IT} = true  # We'll check the actual instance

# For actual PSF instances
supports_3d(psf::P) where P <: AbstractPSF = supports_3d(typeof(psf))
supports_3d(psf::SplinePSF) = psf.z_range !== nothing

"""
    has_z_coordinate(emitter_type::Type{<:AbstractEmitter}) -> Bool

Check if an emitter type has a z-coordinate field.
"""
# Define a single method with a type parameter and conditional logic inside
function has_z_coordinate(::Type{E}) where E <: AbstractEmitter
    hasfield(E, :z)
end

# Method for actual emitter instances
has_z_coordinate(emitter::E) where E <: AbstractEmitter = has_z_coordinate(typeof(emitter))

"""
    _integrate_pixels_core!(
        result::AbstractMatrix{T},
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        func::Function,
        emitter_x::Real,
        emitter_y::Real;
        sampling::Integer=2
    ) where T

Internal function implementing the core pixel integration logic.

# Arguments
- `result`: Pre-allocated array where integration results will be stored
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `func`: Function to evaluate at PSF coordinates
- `emitter_x`, `emitter_y`: Emitter coordinates in same units as pixel edges
- `sampling`: Subpixel sampling density (default: 2)

# Returns
- `result` array filled with integration values
"""
function _integrate_pixels_core!(
    result::AbstractMatrix{T},
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    func::Function,
    emitter_x::Real,
    emitter_y::Real;
    sampling::Integer=2
) where T
    nx = length(pixel_edges_x) - 1
    ny = length(pixel_edges_y) - 1
    
    # Validate dimensions
    size(result) == (ny, nx) || 
        throw(DimensionMismatch("Result array dimensions $(size(result)) do not match expected ($ny, $nx)"))
    
    # Always use multi-threading - has minimal overhead with single thread
    Threads.@threads for ix in 1:nx
        x_start = pixel_edges_x[ix]
        x_end = pixel_edges_x[ix + 1]
        dx = x_end - x_start
        
        for iy in 1:ny
            y_start = pixel_edges_y[iy]
            y_end = pixel_edges_y[iy + 1]
            dy = y_end - y_start
            
            if sampling == 1
                # Just use pixel centers
                x = (x_start + x_end)/2
                y = (y_start + y_end)/2
                result[iy, ix] = func(x - emitter_x, y - emitter_y) * (dx * dy)
            else
                # Integration with subpixel sampling
                dx_sample = dx / sampling
                dy_sample = dy / sampling
                area_factor = dx_sample * dy_sample
                
                # Subsampling points
                xs = range(x_start + dx_sample/2, x_end - dx_sample/2, length=sampling)
                ys = range(y_start + dy_sample/2, y_end - dy_sample/2, length=sampling)
                
                pixel_sum = zero(T)
                for x in xs, y in ys
                    x_psf = x - emitter_x
                    y_psf = y - emitter_y
                    pixel_sum += func(x_psf, y_psf)
                end
                result[iy, ix] = pixel_sum * area_factor
            end
        end
    end
    
    return result
end

"""
    _integrate_pixels_generic!(
        result::AbstractMatrix{T},
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter,
        f::Function;
        sampling::Integer=2
    ) where T

Generic wrapper for PSF integration that handles coordinate systems and validation.
Automatically handles z-coordinates when both PSF and emitter support them.

# Arguments
- `result`: Pre-allocated array where results will be stored
- `psf`: Point spread function to integrate
- `camera`: AbstractCamera defining pixel edges
- `emitter`: Emitter with position information
- `f`: Function to evaluate (e.g., (p,x,y) -> p(x,y) for intensity)
- `sampling`: Subpixel sampling density (default: 2)

# Returns
- `result` array filled with integration values
"""
function _integrate_pixels_generic!(
    result::AbstractMatrix{T},
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter,
    f::Function;
    sampling::Integer=2
) where T
    # Determine if we should use z-coordinate based on PSF and emitter capabilities
    if supports_3d(psf) && has_z_coordinate(emitter)
        # For 3D PSF and emitter with z-coordinate
        eval_func = (x, y) -> f(psf, x, y, emitter.z)
    else
        # For 2D PSF or emitter without z-coordinate
        eval_func = (x, y) -> f(psf, x, y)
    end
    
    # Call core integration function
    return _integrate_pixels_core!(
        result,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        eval_func,
        emitter.x,
        emitter.y;
        sampling=sampling
    )
end

"""
    _integrate_pixels_generic(
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter,
        f::Function,  # Function to integrate
        ::Type{T};   # Return type
        sampling::Integer=2
    ) where T

Internal generic integration routine used by both intensity and amplitude integration.
Creates a new result array and delegates to the in-place version.

# Returns
- Matrix{T} with dimensions [y,x] starting at [1,1] for top-left pixel
"""
function _integrate_pixels_generic(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter,
    f::Function,
    ::Type{T};
    sampling::Integer=2
) where T
    # Get pixel dimensions
    nx = length(camera.pixel_edges_x) - 1
    ny = length(camera.pixel_edges_y) - 1
    
    # Allocate result matrix
    result = Matrix{T}(undef, ny, nx)
    
    # Use new in-place function
    _integrate_pixels_generic!(
        result,
        psf,
        camera,
        emitter,
        f,
        sampling=sampling
    )
    
    return result
end