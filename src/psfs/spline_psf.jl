# src/psfs/spline_psf.jl


"""
    SplinePSF{T<:AbstractFloat, IT<:AbstractInterpolation} <: AbstractPSF

A point spread function (PSF) represented as a B-spline interpolation.

# Fields
- `spline`: The B-spline interpolation object 
- `x_range`: Range of x-coordinates used for uniform grid interpolation
- `y_range`: Range of y-coordinates used for uniform grid interpolation  
- `z_range`: Range of z-coordinates for 3D PSFs, or `nothing` for 2D PSFs
"""
struct SplinePSF{T<:AbstractFloat, IT<:AbstractInterpolation} <: AbstractPSF
    spline::IT
    x_range::StepRangeLen{T}
    y_range::StepRangeLen{T}
    z_range::Union{StepRangeLen{T}, Nothing}
end

# --- Positional argument constructors for AbstractPSF types ---

# Constructor from AbstractPSF with 3D sampling
"""
    SplinePSF(psf::Abstract3DPSF, 
              x_range::AbstractRange, 
              y_range::AbstractRange,
              z_range::AbstractRange;
              order::Integer=3)

Create a 3D SplinePSF by sampling a 3D PSF on a regular grid.
"""
function SplinePSF(psf::Abstract3DPSF, 
                  x_range::AbstractRange,
                  y_range::AbstractRange,
                  z_range::AbstractRange;
                  order::Integer=3)
    # Sample the PSF on the grid
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range), length(z_range))
    for (iz, z) in enumerate(z_range)
        for (ix, x) in enumerate(x_range)
            for (iy, y) in enumerate(y_range)
                psf_stack[iy, ix, iz] = psf(x, y, z)
            end
        end
    end
    
    # Call the array constructor with the sampled data
    return SplinePSF(psf_stack, x_range, y_range, z_range; order=order)
end

# Constructor from AbstractPSF with 2D sampling
"""
    SplinePSF(psf::Abstract2DPSF, 
              x_range::AbstractRange, 
              y_range::AbstractRange;
              order::Integer=3)

Create a 2D SplinePSF by sampling a 2D PSF on a regular grid.
"""
function SplinePSF(psf::Abstract2DPSF, 
                  x_range::AbstractRange,
                  y_range::AbstractRange;
                  order::Integer=3)
    # Sample the PSF on the grid
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range))
    for (ix, x) in enumerate(x_range)
        for (iy, y) in enumerate(y_range)
            psf_stack[iy, ix] = psf(x, y)
        end
    end
    
    # Call the array constructor with the sampled data
    return SplinePSF(psf_stack, x_range, y_range; order=order)
end

# Fallback for any AbstractPSF with 3D sampling
"""
    SplinePSF(psf::AbstractPSF, 
              x_range::AbstractRange, 
              y_range::AbstractRange,
              z_range::AbstractRange;
              order::Integer=3)

Create a 3D SplinePSF by sampling any PSF type on a regular grid.
"""
function SplinePSF(psf::AbstractPSF, 
                  x_range::AbstractRange,
                  y_range::AbstractRange,
                  z_range::AbstractRange;
                  order::Integer=3)
    # Sample the PSF on the grid
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range), length(z_range))
    for (iz, z) in enumerate(z_range)
        for (ix, x) in enumerate(x_range)
            for (iy, y) in enumerate(y_range)
                psf_stack[iy, ix, iz] = psf(x, y, z)
            end
        end
    end
    
    # Call the array constructor with the sampled data
    return SplinePSF(psf_stack, x_range, y_range, z_range; order=order)
end

# Fallback for any AbstractPSF with 2D sampling
"""
    SplinePSF(psf::AbstractPSF, 
              x_range::AbstractRange, 
              y_range::AbstractRange;
              order::Integer=3)

Create a 2D SplinePSF by sampling any PSF type on a regular grid.
"""
function SplinePSF(psf::AbstractPSF, 
                  x_range::AbstractRange,
                  y_range::AbstractRange;
                  order::Integer=3)
    # Sample the PSF on the grid
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range))
    for (ix, x) in enumerate(x_range)
        for (iy, y) in enumerate(y_range)
            psf_stack[iy, ix] = psf(x, y)
        end
    end
    
    # Call the array constructor with the sampled data
    return SplinePSF(psf_stack, x_range, y_range; order=order)
end

# --- Array data constructors ---

# 3D array constructor with positional arguments
"""
    SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
              x_range::AbstractRange,
              y_range::AbstractRange,
              z_range::AbstractRange;
              order::Integer=3)

Construct a 3D SplinePSF from a PSF stack and coordinate ranges.

# Arguments
- `psf_stack`: 3D array containing the PSF data
- `x_range`: Range of uniformly spaced x coordinates
- `y_range`: Range of uniformly spaced y coordinates
- `z_range`: Range of uniformly spaced z coordinates
- `order`: Interpolation order (default: 3)

# Returns
- `SplinePSF`: A 3D spline interpolation of the PSF
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
                  x_range::AbstractRange,
                  y_range::AbstractRange,
                  z_range::AbstractRange;
                  order::Integer=3)
    # Normalize the PSF so it sums to 1
    normalized_psf = psf_stack / sum(psf_stack)
    
    # Convert ranges to StepRangeLen{Float64} for consistent type
    x_range_f64 = convert(StepRangeLen{Float64}, x_range)
    y_range_f64 = convert(StepRangeLen{Float64}, y_range)
    z_range_f64 = convert(StepRangeLen{Float64}, z_range)
    
    # Create the BSpline cubic interpolant
    itp = interpolate(normalized_psf, BSpline(Cubic(Line(OnGrid()))))
    scaled_itp = scale(itp, (y_range_f64, x_range_f64, z_range_f64))
    
    return SplinePSF{Float64, typeof(scaled_itp)}(scaled_itp, x_range_f64, y_range_f64, z_range_f64)
end

# 2D array constructor with positional arguments
"""
    SplinePSF(psf_stack::AbstractArray{<:Real,2}, 
              x_range::AbstractRange,
              y_range::AbstractRange;
              order::Integer=3)

Construct a 2D SplinePSF from a PSF stack and coordinate ranges.

# Arguments
- `psf_stack`: 2D array containing the PSF data
- `x_range`: Range of uniformly spaced x coordinates
- `y_range`: Range of uniformly spaced y coordinates
- `order`: Interpolation order (default: 3)

# Returns
- `SplinePSF`: A 2D spline interpolation of the PSF
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,2}, 
                  x_range::AbstractRange,
                  y_range::AbstractRange;
                  order::Integer=3)
    normalized_psf = psf_stack / sum(psf_stack)
    
    # Convert ranges to StepRangeLen{Float64}
    x_range_f64 = convert(StepRangeLen{Float64}, x_range)
    y_range_f64 = convert(StepRangeLen{Float64}, y_range)
    
    itp = interpolate(normalized_psf, BSpline(Cubic(Line(OnGrid()))))
    scaled_itp = scale(itp, (y_range_f64, x_range_f64))
    
    return SplinePSF{Float64, typeof(scaled_itp)}(scaled_itp, x_range_f64, y_range_f64, nothing)
end

# --- Backward compatibility for keyword argument versions ---

"""
    SplinePSF(psf_stack::AbstractArray{<:Real,3};
              x_range::AbstractRange, 
              y_range::AbstractRange, 
              z_range::AbstractRange,
              order::Integer=3, 
              smoothing::Real=0.0)

Backward compatibility constructor for 3D PSFs.
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,3};
                  x_range::AbstractRange,
                  y_range::AbstractRange,
                  z_range::AbstractRange,
                  order::Integer=3,
                  smoothing::Real=0.0)
    if smoothing > 0
        @warn "Smoothing parameter is ignored; use the positional argument constructor directly."
    end
    return SplinePSF(psf_stack, x_range, y_range, z_range; order=order)
end

"""
    SplinePSF(psf_stack::AbstractArray{<:Real,2};
              x_range::AbstractRange, 
              y_range::AbstractRange,
              order::Integer=3, 
              smoothing::Real=0.0)

Backward compatibility constructor for 2D PSFs.
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,2};
                  x_range::AbstractRange,
                  y_range::AbstractRange,
                  order::Integer=3,
                  smoothing::Real=0.0)
    if smoothing > 0
        @warn "Smoothing parameter is ignored; use the positional argument constructor directly."
    end
    return SplinePSF(psf_stack, x_range, y_range; order=order)
end

# --- Backward compatibility for vector inputs ---
"""
    SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
              x_coords::AbstractVector,
              y_coords::AbstractVector,
              z_coords::AbstractVector;
              order::Integer=3)

Construct a 3D SplinePSF from a PSF stack and coordinate vectors.
This method converts vectors to ranges assuming uniform spacing.
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
                  x_coords::AbstractVector,
                  y_coords::AbstractVector,
                  z_coords::AbstractVector;
                  order::Integer=3)
    @warn "Vector-based constructor is deprecated. Use ranges for uniform grid sampling instead."
    
    # Convert vectors to ranges
    x_range = range(minimum(x_coords), maximum(x_coords), length=length(x_coords))
    y_range = range(minimum(y_coords), maximum(y_coords), length=length(y_coords))
    z_range = range(minimum(z_coords), maximum(z_coords), length=length(z_coords))
    
    return SplinePSF(psf_stack, x_range, y_range, z_range; order=order)
end

"""
    SplinePSF(psf_stack::AbstractArray{<:Real,2};
              x_coords::AbstractVector,
              y_coords::AbstractVector;
              order::Integer=3,
              smoothing::Real=0.0)

Construct a 2D SplinePSF from a PSF stack and coordinate vectors.
This method converts vectors to ranges assuming uniform spacing.

# Deprecated
Use the range-based constructor instead.
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,2};
                   x_coords::AbstractVector,
                   y_coords::AbstractVector;
                   order::Integer=3,
                   smoothing::Real=0.0)
    @warn "Vector-based constructor is deprecated. Use ranges for uniform grid sampling instead."
    
    # Convert vectors to ranges
    x_range = range(minimum(x_coords), maximum(x_coords), length=length(x_coords))
    y_range = range(minimum(y_coords), maximum(y_coords), length=length(y_coords))
    
    return SplinePSF(psf_stack; x_range=x_range, y_range=y_range,
                    order=order, smoothing=smoothing)
end

# --- Convenience constructor from an existing PSF ---


# --- Define evaluation for 3D PSF ---
"""
    (psf::SplinePSF)(x::Real, y::Real, z::Real)

Evaluate the 3D spline PSF at position (x, y, z).

Returns 0 if the position is outside the PSF boundary.
"""
function (psf::SplinePSF)(x::Real, y::Real, z::Real)
    # Check bounds in x and y (and z if applicable)
    if x < first(psf.x_range) || x > last(psf.x_range) ||
       y < first(psf.y_range) || y > last(psf.y_range)
        return zero(Float64)
    end
    if psf.z_range !== nothing && (z < first(psf.z_range) || z > last(psf.z_range))
        return zero(Float64)
    end
    # Evaluate using the scaled interpolant.
    # Note the order: since our data is (y, x, z), we call it as (y, x, z).
    return psf.spline(y, x, z)
end

# --- Define evaluation for 2D PSF ---
"""
    (psf::SplinePSF)(x::Real, y::Real)

Evaluate the 2D spline PSF at position (x, y).

Returns 0 if the position is outside the PSF boundary.
"""
function (psf::SplinePSF)(x::Real, y::Real)
    if x < first(psf.x_range) || x > last(psf.x_range) ||
       y < first(psf.y_range) || y > last(psf.y_range)
        return zero(Float64)
    end
    return psf.spline(y, x)
end

# --- Amplitude functions (unchanged) ---
"""
    amplitude(psf::SplinePSF, x::Real, y::Real, z::Real)

Calculate the complex amplitude of the 3D PSF at position (x, y, z).

Returns the square root of the PSF intensity value.
"""
function amplitude(psf::SplinePSF, x::Real, y::Real, z::Real)
    intensity = psf(x, y, z)
    return sqrt(complex(intensity))
end

"""
    amplitude(psf::SplinePSF, x::Real, y::Real)

Calculate the complex amplitude of the PSF at position (x, y) with z=0.

Returns the square root of the PSF intensity value.
"""
function amplitude(psf::SplinePSF, x::Real, y::Real)
    return amplitude(psf, x, y, 0.0)
end

# Integration methods for pixel calculation

"""
    integrate_pixels(psf::SplinePSF, 
                    camera::AbstractCamera, 
                    emitter::AbstractEmitter;
                    sampling::Integer=2)

Integrate PSF over camera pixels using interpolation.
"""
function integrate_pixels(
    psf::SplinePSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    result = _integrate_pixels_generic(
        psf, camera, emitter, 
        (p, x, y) -> p(x, y, emitter.z),
        Float64; sampling=sampling
    )
    return result ./ sum(result)
end

"""
    integrate_pixels_amplitude(psf::SplinePSF,
                              camera::AbstractCamera,
                              emitter::AbstractEmitter;
                              sampling::Integer=2)

Integrate PSF amplitude (complex) over camera pixels.
"""
function integrate_pixels_amplitude(
    psf::SplinePSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    return _integrate_pixels_generic(
        psf, camera, emitter,
        (p, x, y) -> amplitude(p, x, y, emitter.z),
        Complex{Float64}; sampling=sampling
    )
end

# I/O methods

"""
    save_spline_psf(filename::String, psf::SplinePSF)

Save a SplinePSF to an HDF5 file.
"""
function save_spline_psf(filename::String, psf::SplinePSF)
    h5open(filename, "w") do file
        # Save range information
        file["x_start"] = first(psf.x_range)
        file["x_step"] = step(psf.x_range)
        file["x_length"] = length(psf.x_range)
        
        file["y_start"] = first(psf.y_range)
        file["y_step"] = step(psf.y_range)
        file["y_length"] = length(psf.y_range)
        
        if psf.z_range !== nothing
            file["z_start"] = first(psf.z_range)
            file["z_step"] = step(psf.z_range)
            file["z_length"] = length(psf.z_range)
            
            # Sample the PSF on the grid
            psf_values = Array{Float64}(undef, length(psf.y_range), length(psf.x_range), length(psf.z_range))
            for (iz, z) in enumerate(psf.z_range)
                for (ix, x) in enumerate(psf.x_range)
                    for (iy, y) in enumerate(psf.y_range)
                        psf_values[iy, ix, iz] = psf(x, y, z)
                    end
                end
            end
        else
            # 2D case
            file["z_range"] = "none"
            
            # Sample the PSF on the grid
            psf_values = Array{Float64}(undef, length(psf.y_range), length(psf.x_range))
            for (ix, x) in enumerate(psf.x_range)
                for (iy, y) in enumerate(psf.y_range)
                    psf_values[iy, ix] = psf(x, y)
                end
            end
        end
        
        file["psf_values"] = psf_values
        attrs = file["psf_values"]
        attrs["type"] = "SplinePSF"
        attrs["version"] = "1.1"  # Updated version
    end
end

"""
    load_spline_psf(filename::String)

Load a SplinePSF from an HDF5 file.
"""
function load_spline_psf(filename::String)
    h5open(filename, "r") do file
        # Handle both old (vector-based) and new (range-based) formats
        attrs = attributes(file["psf_values"])
        version = read(attrs["version"])
        
        if version == "1.0"
            # Old format with vectors
            x_coords = read(file["x_knots"])
            y_coords = read(file["y_knots"])
            z_coords = read(file["z_knots"])
            psf_values = read(file["psf_values"])
            
            # Convert to ranges
            x_range = range(minimum(x_coords), maximum(x_coords), length=length(x_coords))
            y_range = range(minimum(y_coords), maximum(y_coords), length=length(y_coords))
            z_range = range(minimum(z_coords), maximum(z_coords), length=length(z_coords))
            
            return SplinePSF(psf_values; 
                            x_range=x_range, 
                            y_range=y_range, 
                            z_range=z_range)
        else
            # New format with ranges
            x_start = read(file["x_start"])
            x_step = read(file["x_step"])
            x_length = read(file["x_length"])
            
            y_start = read(file["y_start"])
            y_step = read(file["y_step"])
            y_length = read(file["y_length"])
            
            psf_values = read(file["psf_values"])
            
            # Create ranges from saved parameters
            x_range = range(x_start, step=x_step, length=x_length)
            y_range = range(y_start, step=y_step, length=y_length)
            
            # Check if it's a 3D PSF
            if haskey(file, "z_start")
                z_start = read(file["z_start"])
                z_step = read(file["z_step"])
                z_length = read(file["z_length"])
                z_range = range(z_start, step=z_step, length=z_length)
                
                return SplinePSF(psf_values; 
                                x_range=x_range, 
                                y_range=y_range, 
                                z_range=z_range)
            else
                # 2D PSF
                return SplinePSF(psf_values; 
                                x_range=x_range, 
                                y_range=y_range)
            end
        end
    end
end

# Pretty printing
function Base.show(io::IO, psf::SplinePSF)
    nx = length(psf.x_range)
    ny = length(psf.y_range)
    
    x_size = last(psf.x_range) - first(psf.x_range)
    
    if psf.z_range !== nothing
        nz = length(psf.z_range)
        z_size = last(psf.z_range) - first(psf.z_range)
        print(io, "SplinePSF($(ny)×$(nx)×$(nz) grid, $(round(x_size, digits=2))×$(round(z_size, digits=2))μm)")
    else
        print(io, "SplinePSF($(ny)×$(nx) grid, $(round(x_size, digits=2))μm)")
    end
end
