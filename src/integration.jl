# src/integration.jl

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
    
    # Calculate pixel dimensions and area factors
    dx = camera.pixel_edges_x[2] - camera.pixel_edges_x[1]
    dy = camera.pixel_edges_y[2] - camera.pixel_edges_y[1]
    pixel_area = dx * dy
    
    result = Matrix{T}(undef, ny, nx)
    
    # For each pixel
    Threads.@threads for ix in 1:nx
        for iy in 1:ny
            # Get physical coordinates of pixel edges
            x_start = camera.pixel_edges_x[ix]
            x_end = camera.pixel_edges_x[ix + 1]
            y_start = camera.pixel_edges_y[iy]
            y_end = camera.pixel_edges_y[iy + 1]
            
            if sampling == 1
                # Just use pixel centers
                x = (x_start + x_end)/2
                y = (y_start + y_end)/2
                result[iy, ix] = f(psf, x - emitter.x, y - emitter.y) * pixel_area
            else
                # Integration step size
                dx_sample = dx / sampling
                dy_sample = dy / sampling
                area_factor = dx_sample * dy_sample
                
                # Subsampling points
                xs = range(x_start + dx_sample/2, x_end - dx_sample/2, sampling)
                ys = range(y_start + dy_sample/2, y_end - dy_sample/2, sampling)
                
                pixel_sum = zero(T)
                for x in xs, y in ys
                    x_psf = x - emitter.x
                    y_psf = y - emitter.y
                    pixel_sum += f(psf, x_psf, y_psf)
                end
                result[iy, ix] = pixel_sum * area_factor
            end
        end
    end
    
    return result
end

"""
    integrate_pixels(psf::AbstractPSF, camera::AbstractCamera, emitter::AbstractEmitter; sampling::Integer=2)

Integrate PSF intensity over camera pixels.

For each pixel in the camera, numerically integrates the PSF intensity using the specified 
sampling density. Physical coordinates are relative to camera with (0,0) at top-left corner.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `camera::AbstractCamera`: Camera geometry defining pixel edges in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `sampling::Integer=2`: Number of samples per pixel in each dimension

# Returns
- `Matrix{T}`: Integrated intensities where T matches emitter.photons type
- Array is indexed as [y,x] with [1,1] at top-left pixel
- Values are normalized to sum to 1

# Examples
```julia
camera = IdealCamera(0:0.1:2.0, 0:0.1:2.0)  # 20x20 camera, 100nm pixels
emitter = Emitter2D(1.0, 1.0, 1000.0)       # Emitter at (1μm, 1μm)
psf = Gaussian2D(0.15)                       # σ = 150nm
pixels = integrate_pixels(psf, camera, emitter)
```

See also: [`integrate_pixels_amplitude`](@ref), [`AbstractPSF`](@ref)
"""
function integrate_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    T = typeof(emitter.photons)
    result = _integrate_pixels_generic(psf, camera, emitter, (p,x,y) -> p(x,y), T; sampling=sampling)
    return result ./ sum(result)
end

"""
    integrate_pixels_amplitude(psf::AbstractPSF, camera::AbstractCamera, emitter::AbstractEmitter; sampling::Integer=2)

Integrate PSF complex amplitude over camera pixels.

For each pixel in the camera, numerically integrates the PSF amplitude using the specified
sampling density. Unlike intensity integration, returns unnormalized complex amplitudes
which can be used for coherent calculations.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `camera::AbstractCamera`: Camera geometry defining pixel edges in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `sampling::Integer=2`: Number of samples per pixel in each dimension

# Returns
- `Matrix{Complex{T}}`: Integrated complex amplitudes where T matches emitter.photons type
- Array is indexed as [y,x] with [1,1] at top-left pixel
- Values are not normalized to preserve complex amplitude relationships

# Notes
- For coherent calculations, use this function instead of `integrate_pixels`
- Return type is complex to support PSFs with phase information
- To get intensity from amplitude: `abs2.(integrate_pixels_amplitude(...))`

# Examples
```julia
camera = IdealCamera(0:0.1:2.0, 0:0.1:2.0)  # 20x20 camera, 100nm pixels
emitter = Emitter2D(1.0, 1.0, 1000.0)       # Emitter at (1μm, 1μm)
psf = Gaussian2D(0.15)                       # σ = 150nm
amplitudes = integrate_pixels_amplitude(psf, camera, emitter)
intensities = abs2.(amplitudes)              # Convert to intensity if needed
```

See also: [`integrate_pixels`](@ref), [`AbstractPSF`](@ref)
"""
function integrate_pixels_amplitude(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    T = Complex{typeof(emitter.photons)}
    return _integrate_pixels_generic(psf, camera, emitter, amplitude, T; sampling=sampling)
end