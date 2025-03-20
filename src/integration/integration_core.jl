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
    get_pixel_indices(
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter,
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}
    ) -> Tuple{UnitRange{Int},UnitRange{Int}}

Get the pixel indices that overlap with the support region.

# Arguments
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `emitter`: Emitter with position information
- `support`: Region of interest, either a radius or explicit bounds (x_min, x_max, y_min, y_max)

# Returns
- Tuple of (i_range, j_range) with pixel indices that cover the support region
"""
function get_pixel_indices(
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter,
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}
)
    if support isa Real
        # Convert radius to rectangular region
        if isinf(support)
            # Use full grid
            return 1:length(pixel_edges_x)-1, 1:length(pixel_edges_y)-1
        else
            # Define region around emitter
            x_min = emitter.x - support
            x_max = emitter.x + support
            y_min = emitter.y - support
            y_max = emitter.y + support
        end
    else
        # Explicit region provided as tuple
        x_min, x_max, y_min, y_max = support
    end
    
    # Find pixel indices that overlap with region
    i_min = searchsortedlast(pixel_edges_x, x_min)
    i_max = searchsortedfirst(pixel_edges_x, x_max)
    j_min = searchsortedlast(pixel_edges_y, y_min)
    j_max = searchsortedfirst(pixel_edges_y, y_max)
    
    # Ensure at least one pixel in each dimension
    i_min = max(1, i_min)
    i_max = min(length(pixel_edges_x)-1, max(i_min, i_max))
    j_min = max(1, j_min)
    j_max = min(length(pixel_edges_y)-1, max(j_min, j_max))
    
    return i_min:i_max, j_min:j_max
end

"""
    _integrate_pixels_core!(
        result::AbstractArray,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        func::Function,
        emitter_x::Real,
        emitter_y::Real;
        sampling::Integer=2,
        threaded::Bool=true
    )

Internal function implementing the core pixel integration logic.
Supports both scalar and vector/tensor return types from the function.

# Arguments
- `result`: Pre-allocated array where integration results will be stored
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `func`: Function to evaluate at PSF coordinates - can return scalar, vector, or tensor
- `emitter_x`, `emitter_y`: Emitter coordinates in same units as pixel edges
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)
  Set to false when using with automatic differentiation frameworks

# Returns
- `result` array filled with integration values
"""
function _integrate_pixels_core!(
    result::AbstractArray,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    func::Function,
    emitter_x::Real,
    emitter_y::Real;
    sampling::Integer=2,
    threaded::Bool=true
)
    nx = length(pixel_edges_x) - 1
    ny = length(pixel_edges_y) - 1
    
    # Get basic result dimensions and validate
    if result isa AbstractMatrix
        # For 2D result (scalar return case)
        size(result) == (ny, nx) || 
            throw(DimensionMismatch("Result array dimensions $(size(result)) do not match expected ($ny, $nx)"))
    elseif result isa AbstractArray{<:Any,3}
        # For 3D result (vector return case)
        size(result)[1:2] == (ny, nx) || 
            throw(DimensionMismatch("First two dimensions of result $(size(result)[1:2]) do not match expected ($ny, $nx)"))
    else
        throw(ArgumentError("Unsupported result array type: $(typeof(result))"))
    end
    
    # Define a nested function to process a single pixel
    function process_pixel(ix, iy)
        x_start = pixel_edges_x[ix]
        x_end = pixel_edges_x[ix + 1]
        dx = x_end - x_start
        
        y_start = pixel_edges_y[iy]
        y_end = pixel_edges_y[iy + 1]
        dy = y_end - y_start
        
        # Special case for single sampling point
        if sampling == 1
            # Just use pixel centers
            x = (x_start + x_end)/2
            y = (y_start + y_end)/2
            
            # Calculate the function value once
            func_val = func(x - emitter_x, y - emitter_y)
            area = dx * dy
            
            # Handle different return types
            if func_val isa Number
                # Scalar case
                result[iy, ix] = func_val * area
            else
                # Vector or tensor case - use broadcasting to scale all elements
                result[iy, ix, :] .= func_val .* area
            end
        else
            # Integration with subpixel sampling
            dx_sample = dx / sampling
            dy_sample = dy / sampling
            area_factor = dx_sample * dy_sample
            
            # Subsampling points
            xs = range(x_start + dx_sample/2, x_end - dx_sample/2, length=sampling)
            ys = range(y_start + dy_sample/2, y_end - dy_sample/2, length=sampling)
            
            # Get first value to determine return type and initialize sum properly
            x_psf = xs[1] - emitter_x
            y_psf = ys[1] - emitter_y
            first_val = func(x_psf, y_psf)
            
            if first_val isa Number
                # Scalar case
                pixel_sum = zero(typeof(first_val))
                for x in xs, y in ys
                    pixel_sum += func(x - emitter_x, y - emitter_y)
                end
                result[iy, ix] = pixel_sum * area_factor
            else
                # Vector case - initialize sum with correct dimensions
                pixel_sum = zeros(eltype(first_val), length(first_val))
                for x in xs, y in ys
                    pixel_sum .+= func(x - emitter_x, y - emitter_y)
                end
                # Store each component scaled by area
                result[iy, ix, :] .= pixel_sum .* area_factor
            end
        end
    end
    
    # Apply with or without threading
    if threaded
        Threads.@threads for ix in 1:nx
            for iy in 1:ny
                process_pixel(ix, iy)
            end
        end
    else
        for ix in 1:nx
            for iy in 1:ny
                process_pixel(ix, iy)
            end
        end
    end
    
    return result
end

"""
    _integrate_pixels_generic!(
        result::AbstractArray,
        psf::AbstractPSF,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter,
        f::Function;
        sampling::Integer=2,
        threaded::Bool=true
    )

Generic wrapper for PSF integration that handles coordinate systems and validation.
Automatically handles z-coordinates when both PSF and emitter support them.
Supports functions that return either scalars or vectors/tensors.

# Arguments
- `result`: Pre-allocated array where results will be stored
- `psf`: Point spread function to integrate
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `emitter`: Emitter with position information
- `f`: Function to evaluate (e.g., (p,x,y) -> p(x,y) for intensity)
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- `result` array filled with integration values
"""
function _integrate_pixels_generic!(
    result::AbstractArray,
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter,
    f::Function;
    sampling::Integer=2,
    threaded::Bool=true
)
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
        pixel_edges_x,
        pixel_edges_y,
        eval_func,
        emitter.x,
        emitter.y;
        sampling=sampling,
        threaded=threaded
    )
end

"""
    _integrate_pixels_generic(
        psf::AbstractPSF,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter,
        f::Function,
        ::Type{T},
        output_dims=();
        sampling::Integer=2,
        threaded::Bool=true
    ) where T

Internal generic integration routine used by both intensity and amplitude integration.
Creates a new result array and delegates to the in-place version.
Supports both scalar and vector/tensor return types.

# Arguments
- `output_dims`: Additional dimensions for the result array beyond the basic [y,x] dimensions.
  For vector returns (e.g., [Ex, Ey]), set output_dims=(2,) for a [y,x,2] result array.
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- Array with dimensions [y,x] or [y,x,output_dims...] depending on the function's return type
"""
function _integrate_pixels_generic(
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter,
    f::Function,
    ::Type{T},
    output_dims=();
    sampling::Integer=2,
    threaded::Bool=true
) where T
    # Get pixel dimensions
    nx = length(pixel_edges_x) - 1
    ny = length(pixel_edges_y) - 1
    
    # Allocate result array with appropriate dimensions
    if isempty(output_dims)
        # Scalar case - 2D result
        result = Array{T}(undef, ny, nx)
    else
        # Vector/tensor case - 3D or higher result
        result = Array{T}(undef, ny, nx, output_dims...)
    end
    
    # Use in-place function
    _integrate_pixels_generic!(
        result,
        psf,
        pixel_edges_x,
        pixel_edges_y,
        emitter,
        f,
        sampling=sampling,
        threaded=threaded
    )
    
    return result
end