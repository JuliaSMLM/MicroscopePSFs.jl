# src/psfs/spline_psf.jl

"""
    SplinePSF{T<:AbstractFloat, IT<:AbstractInterpolation} <: AbstractPSF

A point spread function (PSF) represented as a B-spline interpolation.

# Fields
- `spline`: The B-spline interpolation object 
- `x_range`: Range of x-coordinates used for uniform grid interpolation
- `y_range`: Range of y-coordinates used for uniform grid interpolation  
- `z_range`: Range of z-coordinates for 3D PSFs, or `nothing` for 2D PSFs

# Notes
- Coordinates and ranges are in physical units (typically microns)
- PSF values are preserved from the original PSF that was sampled
- Internally uses cubic B-spline interpolation (order=3) by default
"""
struct SplinePSF{T<:AbstractFloat, IT<:AbstractInterpolation} <: AbstractPSF
    spline::IT
    x_range::StepRangeLen{T}
    y_range::StepRangeLen{T}
    z_range::Union{StepRangeLen{T}, Nothing}
end

# --- Helper functions for PSF sampling ---

"""
    _sample_psf_3d(psf::AbstractPSF, x_range::AbstractRange, y_range::AbstractRange, z_range::AbstractRange)

Internal helper function to sample a 3D PSF on a regular grid.
Returns an array with dimensions [y, x, z] containing raw PSF values.
"""
function _sample_psf_3d(psf::AbstractPSF, x_range::AbstractRange, y_range::AbstractRange, z_range::AbstractRange)
    # Sample the PSF on the grid
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range), length(z_range))
    for (iz, z) in enumerate(z_range)
        for (ix, x) in enumerate(x_range)
            for (iy, y) in enumerate(y_range)
                psf_stack[iy, ix, iz] = psf(x, y, z)
            end
        end
    end
    
    # Return the raw samples without normalization
    return psf_stack
end

"""
    _sample_psf_2d(psf::AbstractPSF, x_range::AbstractRange, y_range::AbstractRange)

Internal helper function to sample a 2D PSF on a regular grid.
Returns an array with dimensions [y, x] containing raw PSF values.
"""
function _sample_psf_2d(psf::AbstractPSF, x_range::AbstractRange, y_range::AbstractRange)
    # Sample the PSF on the grid
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range))
    for (ix, x) in enumerate(x_range)
        for (iy, y) in enumerate(y_range)
            psf_stack[iy, ix] = psf(x, y)
        end
    end
    
    # Return the raw samples without normalization
    return psf_stack
end

# --- Core constructors for 3D PSFs ---

"""
    SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
              x_range::AbstractRange,
              y_range::AbstractRange,
              z_range::AbstractRange;
              order::Integer=3)

Construct a 3D SplinePSF from a PSF stack and coordinate ranges.

# Arguments
- `psf_stack`: 3D array containing the PSF data with dimensions [y, x, z]
- `x_range`: Range of uniformly spaced x coordinates in microns
- `y_range`: Range of uniformly spaced y coordinates in microns
- `z_range`: Range of uniformly spaced z coordinates in microns
- `order`: Interpolation order (default: 3 for cubic B-splines)

# Returns
- `SplinePSF`: A 3D spline interpolation of the PSF
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
                  x_range::AbstractRange,
                  y_range::AbstractRange,
                  z_range::AbstractRange;
                  order::Integer=3)
    # Convert ranges to StepRangeLen{Float64} for consistent type
    x_range_f64 = convert(StepRangeLen{Float64}, x_range)
    y_range_f64 = convert(StepRangeLen{Float64}, y_range)
    z_range_f64 = convert(StepRangeLen{Float64}, z_range)
    
    # Create the BSpline interpolant with the specified order
    if order == 3  # Cubic interpolation (most common case)
        itp = interpolate(psf_stack, BSpline(Cubic(Line(OnGrid()))))
    elseif order == 1  # Linear interpolation
        itp = interpolate(psf_stack, BSpline(Linear()))
    elseif order == 0  # Constant interpolation
        itp = interpolate(psf_stack, BSpline(Constant()))
    else
        throw(ArgumentError("Interpolation order $order not supported. Use 0 (constant), 1 (linear), or 3 (cubic)"))
    end
    
    # Scale the interpolant to the provided coordinate ranges
    # The scale function expects individual ranges, not a tuple
    scaled_itp = scale(itp, y_range_f64, x_range_f64, z_range_f64)
    
    return SplinePSF{Float64, typeof(scaled_itp)}(scaled_itp, x_range_f64, y_range_f64, z_range_f64)
end

# --- Core constructors for 2D PSFs ---

"""
    SplinePSF(psf_stack::AbstractArray{<:Real,2}, 
              x_range::AbstractRange,
              y_range::AbstractRange;
              order::Integer=3)

Construct a 2D SplinePSF from a PSF stack and coordinate ranges.

# Arguments
- `psf_stack`: 2D array containing the PSF data with dimensions [y, x]
- `x_range`: Range of uniformly spaced x coordinates in microns
- `y_range`: Range of uniformly spaced y coordinates in microns
- `order`: Interpolation order (default: 3 for cubic B-splines)

# Returns
- `SplinePSF`: A 2D spline interpolation of the PSF
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,2}, 
                  x_range::AbstractRange,
                  y_range::AbstractRange;
                  order::Integer=3)
    # Convert ranges to StepRangeLen{Float64}
    x_range_f64 = convert(StepRangeLen{Float64}, x_range)
    y_range_f64 = convert(StepRangeLen{Float64}, y_range)
    
    # Create the BSpline interpolant with the specified order
    if order == 3  # Cubic interpolation (most common case)
        itp = interpolate(psf_stack, BSpline(Cubic(Line(OnGrid()))))
    elseif order == 1  # Linear interpolation
        itp = interpolate(psf_stack, BSpline(Linear()))
    elseif order == 0  # Constant interpolation
        itp = interpolate(psf_stack, BSpline(Constant()))
    else
        throw(ArgumentError("Interpolation order $order not supported. Use 0 (constant), 1 (linear), or 3 (cubic)"))
    end
    
    # Scale the interpolant to the provided coordinate ranges
    # The scale function expects individual ranges, not a tuple
    scaled_itp = scale(itp, y_range_f64, x_range_f64)
    
    return SplinePSF{Float64, typeof(scaled_itp)}(scaled_itp, x_range_f64, y_range_f64, nothing)
end

# --- Constructors for building from existing PSF types ---

"""
    SplinePSF(psf::AbstractPSF, 
              x_range::AbstractRange, 
              y_range::AbstractRange,
              z_range::AbstractRange;
              order::Integer=3)

Create a 3D SplinePSF by sampling an existing PSF on a regular grid.

# Arguments
- `psf`: Source PSF to sample
- `x_range`: Range of x-coordinates to sample in microns
- `y_range`: Range of y-coordinates to sample in microns
- `z_range`: Range of z-coordinates to sample in microns
- `order`: Interpolation order (default: 3 for cubic)

# Returns
- `SplinePSF`: A spline interpolation of the sampled PSF

# Example
```julia
# Create a spline from a Scalar3DPSF with 1µm spacing over 4µm range
scalar_psf = Scalar3DPSF(1.4, 0.532, 1.518)  # NA=1.4, λ=532nm, n=1.518
x_range = y_range = range(-2.0, 2.0, length=41)
z_range = range(-1.0, 1.0, length=21)
spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)
```
"""
function SplinePSF(psf::AbstractPSF, 
                  x_range::AbstractRange, 
                  y_range::AbstractRange,
                  z_range::AbstractRange;
                  order::Integer=3)
    # Sample the PSF on the grid
    psf_stack = _sample_psf_3d(psf, x_range, y_range, z_range)
    
    # Create the SplinePSF using the array constructor
    return SplinePSF(psf_stack, x_range, y_range, z_range; order=order)
end

"""
    SplinePSF(psf::AbstractPSF, 
              x_range::AbstractRange, 
              y_range::AbstractRange;
              order::Integer=3)

Create a 2D SplinePSF by sampling a PSF on a 2D grid (at z=0).

# Arguments
- `psf`: Source PSF to sample
- `x_range`: Range of x-coordinates to sample in microns
- `y_range`: Range of y-coordinates to sample in microns
- `order`: Interpolation order (default: 3 for cubic)

# Returns
- `SplinePSF`: A 2D spline interpolation of the sampled PSF

# Note
- For 3D PSFs, this samples at z=0 only
"""
function SplinePSF(psf::AbstractPSF, 
                  x_range::AbstractRange, 
                  y_range::AbstractRange;
                  order::Integer=3)
    # Sample the PSF on the grid (at z=0 for 3D PSFs)
    psf_stack = _sample_psf_2d(psf, x_range, y_range)
    
    # Create the SplinePSF using the array constructor
    return SplinePSF(psf_stack, x_range, y_range; order=order)
end

# --- Convenience constructor with auto-ranges ---

"""
    SplinePSF(psf::AbstractPSF; 
              lateral_range::Float64=2.0,
              axial_range::Float64=1.0,
              lateral_step::Float64=0.05,
              axial_step::Float64=0.1,
              order::Integer=3)

Create a SplinePSF with automatically calculated sampling ranges.

# Arguments
- `psf`: Source PSF to sample
- `lateral_range`: Half-width of lateral (xy) sampling range in microns
- `axial_range`: Half-width of axial (z) sampling range in microns
- `lateral_step`: Step size in microns for lateral (xy) sampling (default: 0.05µm = 50nm)
- `axial_step`: Step size in microns for axial (z) sampling (default: 0.1µm = 100nm)
- `order`: Interpolation order (default: 3 for cubic)

# Returns
- `SplinePSF`: Either a 2D or 3D spline PSF depending on input PSF type

# Example
```julia
# Create a spline PSF with default sampling parameters
scalar_psf = Scalar3DPSF(1.4, 0.532, 1.518)
spline_psf = SplinePSF(scalar_psf)  # Uses default ranges and step sizes
```
"""
function SplinePSF(psf::AbstractPSF; 
                  lateral_range::Float64=2.0,
                  axial_range::Float64=1.0,
                  lateral_step::Float64=0.05,
                  axial_step::Float64=0.1,
                  order::Integer=3)
    # Calculate ranges with specified step size
    x_range = y_range = range(-lateral_range, lateral_range, step=lateral_step)
    
    # Check if this is a 3D PSF by trying to call it with a z parameter
    is_3d = false
    try
        psf(0.0, 0.0, 0.0)
        is_3d = true
    catch e
        if e isa MethodError
            is_3d = false
        else
            rethrow(e)
        end
    end
    
    # Create either a 2D or 3D SplinePSF based on the input PSF
    if is_3d
        z_range = range(-axial_range, axial_range, step=axial_step)
        return SplinePSF(psf, x_range, y_range, z_range; order=order)
    else
        return SplinePSF(psf, x_range, y_range; order=order)
    end
end

# --- Define evaluation for 3D PSF ---

"""
    (psf::SplinePSF)(x::Real, y::Real, z::Real)

Evaluate the 3D spline PSF at position (x, y, z).

# Arguments
- `x, y, z`: Coordinates in microns relative to PSF center

# Returns
- PSF intensity at the specified position

# Notes
- Returns 0 if the position is outside the PSF boundary
- Preserves the normalization of the original PSF data
"""
function (psf::SplinePSF)(x::Real, y::Real, z::Real)
    # Check that this is a 3D PSF
    if psf.z_range === nothing
        throw(ArgumentError("Cannot evaluate a 2D SplinePSF with 3D coordinates. Use psf(x, y) instead."))
    end
    
    # Check bounds in x, y, and z
    if x < first(psf.x_range) || x > last(psf.x_range) ||
       y < first(psf.y_range) || y > last(psf.y_range) ||
       z < first(psf.z_range) || z > last(psf.z_range)
        return zero(Float64)
    end
    
    # Evaluate using the scaled interpolant.
    # Note the order: since our data is (y, x, z), we call it as (y, x, z).
    return psf.spline(y, x, z)
end

# --- Define evaluation for 2D PSF ---

"""
    (psf::SplinePSF)(x::Real, y::Real)

Evaluate the spline PSF at position (x, y).

# Arguments
- `x, y`: Coordinates in microns relative to PSF center

# Returns
- PSF intensity at the specified position

# Notes
- For 3D PSFs, evaluates at z = 0 if in range, otherwise returns 0
- Returns 0 if the position is outside the PSF boundary
- Preserves the normalization of the original PSF data
"""
function (psf::SplinePSF)(x::Real, y::Real)
    # Check bounds in x and y
    if x < first(psf.x_range) || x > last(psf.x_range) ||
       y < first(psf.y_range) || y > last(psf.y_range)
        return zero(Float64)
    end
    
    # For 2D PSFs
    if psf.z_range === nothing
        return psf.spline(y, x)
    else
        # For 3D PSFs, check if z=0 is in range
        if first(psf.z_range) <= 0.0 && last(psf.z_range) >= 0.0
            return psf.spline(y, x, 0.0)
        else
            return zero(Float64)
        end
    end
end

# --- Amplitude functions ---

"""
    amplitude(psf::SplinePSF, x::Real, y::Real, z::Real)

Calculate the complex amplitude of the 3D PSF at position (x, y, z).

# Arguments
- `psf`: SplinePSF instance
- `x, y, z`: Coordinates in microns relative to PSF center

# Returns
- Complex amplitude = sqrt(intensity) with zero phase

# Notes
- Returns sqrt(intensity) as Complex to match the interface of other PSFs
- The SplinePSF does not model phase information
"""
function amplitude(psf::SplinePSF, x::Real, y::Real, z::Real)
    intensity = psf(x, y, z)
    return sqrt(complex(intensity))
end

"""
    amplitude(psf::SplinePSF, x::Real, y::Real)

Calculate the complex amplitude of the PSF at position (x, y) with z=0.

# Arguments
- `psf`: SplinePSF instance
- `x, y`: Coordinates in microns relative to PSF center

# Returns
- Complex amplitude = sqrt(intensity) with zero phase
"""
function amplitude(psf::SplinePSF, x::Real, y::Real)
    intensity = psf(x, y)
    return sqrt(complex(intensity))
end

# --- Integration methods for pixel calculation ---

"""
    integrate_pixels(psf::SplinePSF, 
                    camera::AbstractCamera, 
                    emitter::AbstractEmitter;
                    sampling::Integer=2)

Integrate PSF over camera pixels using interpolation.

# Arguments
- `psf`: SplinePSF instance
- `camera`: Camera geometry
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density for integration accuracy

# Returns
- Array of integrated PSF intensities with dimensions [ny, nx]
- Values represent actual photon counts based on emitter's photon value
"""
function integrate_pixels(
    psf::SplinePSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    # For 3D PSFs using z from the emitter
    if psf.z_range !== nothing
        if !hasfield(typeof(emitter), :z)
            throw(ArgumentError("3D SplinePSF requires an emitter with a z-coordinate"))
        end
        
        result = _integrate_pixels_generic(
            psf, camera, emitter, 
            (p, x, y) -> p(x, y, emitter.z),
            Float64; sampling=sampling
        )
    else
        # For 2D PSFs
        result = _integrate_pixels_generic(
            psf, camera, emitter, 
            (p, x, y) -> p(x, y),
            Float64; sampling=sampling
        )
    end
    
    # Multiply by photon count to preserve physical meaning
    return result .* emitter.photons
end

"""
    integrate_pixels_amplitude(psf::SplinePSF,
                             camera::AbstractCamera,
                             emitter::AbstractEmitter;
                             sampling::Integer=2)

Integrate PSF amplitude (complex) over camera pixels.

# Arguments
- `psf`: SplinePSF instance
- `camera`: Camera geometry
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density for integration accuracy

# Returns
- Array of integrated PSF complex amplitudes with dimensions [ny, nx]
"""
function integrate_pixels_amplitude(
    psf::SplinePSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    # For 3D PSFs using z from the emitter
    if psf.z_range !== nothing
        if !hasfield(typeof(emitter), :z)
            throw(ArgumentError("3D SplinePSF requires an emitter with a z-coordinate"))
        end
        
        return _integrate_pixels_generic(
            psf, camera, emitter, 
            (p, x, y) -> amplitude(p, x, y, emitter.z),
            Complex{Float64}; sampling=sampling
        )
    else
        # For 2D PSFs
        return _integrate_pixels_generic(
            psf, camera, emitter, 
            (p, x, y) -> amplitude(p, x, y),
            Complex{Float64}; sampling=sampling
        )
    end
end

# --- I/O methods ---

"""
    save_spline_psf(filename::String, psf::SplinePSF)

Save a SplinePSF to an HDF5 file.

# Arguments
- `filename`: Path where the PSF will be saved
- `psf`: SplinePSF to save

# Notes
- Saves the coordinates and PSF values to reconstruct the spline later
- File can be loaded with load_spline_psf
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
        
        # Check if 3D
        if psf.z_range !== nothing
            file["z_start"] = first(psf.z_range)
            file["z_step"] = step(psf.z_range)
            file["z_length"] = length(psf.z_range)
            file["dimensions"] = "3D"
            
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
            file["dimensions"] = "2D"
            
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
        attrs["version"] = "2.0"  # Updated version
    end
end

"""
    load_spline_psf(filename::String)

Load a SplinePSF from an HDF5 file.

# Arguments
- `filename`: Path to the saved PSF file

# Returns
- SplinePSF reconstructed from the file data

# Notes
- Handles files created with save_spline_psf
"""
function load_spline_psf(filename::String)
    h5open(filename, "r") do file
        # Read range information
        x_start = read(file["x_start"])
        x_step = read(file["x_step"])
        x_length = read(file["x_length"])
        
        y_start = read(file["y_start"])
        y_step = read(file["y_step"])
        y_length = read(file["y_length"])
        
        # Create ranges from saved parameters
        x_range = range(x_start, step=x_step, length=x_length)
        y_range = range(y_start, step=y_step, length=y_length)
        
        # Read PSF values
        psf_values = read(file["psf_values"])
        
        # Check if it's a 3D PSF
        dimensions = read(file["dimensions"])
        if dimensions == "3D"
            z_start = read(file["z_start"])
            z_step = read(file["z_step"])
            z_length = read(file["z_length"])
            z_range = range(z_start, step=z_step, length=z_length)
            
            return SplinePSF(psf_values, x_range, y_range, z_range)
        else
            # 2D PSF
            return SplinePSF(psf_values, x_range, y_range)
        end
    end
end

# --- Pretty printing ---

"""
    Base.show(io::IO, psf::SplinePSF)

Custom display method for SplinePSF instances.
"""
function Base.show(io::IO, psf::SplinePSF)
    nx = length(psf.x_range)
    ny = length(psf.y_range)
    
    x_size = last(psf.x_range) - first(psf.x_range)
    y_size = last(psf.y_range) - first(psf.y_range)
    
    if psf.z_range !== nothing
        nz = length(psf.z_range)
        z_size = last(psf.z_range) - first(psf.z_range)
        print(io, "SplinePSF($(ny)×$(nx)×$(nz) grid, $(round(x_size, digits=2))×$(round(y_size, digits=2))×$(round(z_size, digits=2))μm)")
    else
        print(io, "SplinePSF($(ny)×$(nx) grid, $(round(x_size, digits=2))×$(round(y_size, digits=2))μm)")
    end
end