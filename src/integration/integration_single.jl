# src/integration/integration_single.jl

"""
    integrate_pixels!(
        result::AbstractMatrix{T},
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2
    ) where T <: Real

Integrate PSF intensity over camera pixels, storing the result in a pre-allocated matrix.
Automatically uses z-coordinate if both PSF and emitter support it.

# Arguments
- `result`: Pre-allocated array where results will be stored
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)

# Returns
- The `result` array, now filled with integrated intensities scaled by emitter.photons

# Notes
- Array is indexed as [y,x] with [1,1] at top-left pixel
"""
function integrate_pixels!(
    result::AbstractMatrix{T},
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
) where T <: Real
    
    # Perform integration
    _integrate_pixels_generic!(
        result,
        psf,
        camera,
        emitter,
        (p, x, y) -> p(x, y),  # For 2D PSF
        sampling=sampling
    )
    
    # Scale by photon count
    result .*= emitter.photons
    
    return result
end

"""
    integrate_pixels_amplitude!(
        result::AbstractMatrix{Complex{T}},
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2
    ) where T <: Real

Integrate PSF complex amplitude over camera pixels, storing the result in a pre-allocated matrix.
Automatically uses z-coordinate if both PSF and emitter support it.

# Arguments
- `result`: Pre-allocated complex array where results will be stored
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)

# Returns
- The `result` array, now filled with integrated complex amplitudes

# Notes
- Array is indexed as [y,x] with [1,1] at top-left pixel
"""
function integrate_pixels_amplitude!(
    result::AbstractMatrix{Complex{T}},
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
) where T <: Real
    
    # Perform integration
    _integrate_pixels_generic!(
        result,
        psf,
        camera,
        emitter,
        amplitude,  # Function that automatically handles 2D or 3D based on context
        sampling=sampling
    )
    
    return result
end

"""
    integrate_pixels(psf::AbstractPSF, camera::AbstractCamera, emitter::AbstractEmitter; sampling::Integer=2)

Integrate PSF intensity over camera pixels.

For each pixel in the camera, numerically integrates the PSF intensity using the specified 
sampling density. Physical coordinates are relative to camera with (0,0) at top-left corner.
Automatically handles z-coordinate if both PSF and emitter support it.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `camera::AbstractCamera`: Camera geometry defining pixel edges in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `sampling::Integer=2`: Number of samples per pixel in each dimension

# Returns
- `Matrix{T}`: Integrated intensities where T matches emitter.photons type
- Array is indexed as [y,x] with [1,1] at top-left pixel
- Values represent actual photon counts in each pixel based on emitter's photon value

# Examples
```julia
camera = IdealCamera(0:0.1:2.0, 0:0.1:2.0)  # 20x20 camera, 100nm pixels
emitter = Emitter2D(1.0, 1.0, 1000.0)       # Emitter at (1μm, 1μm) with 1000 photons
psf = Gaussian2D(0.15)                       # σ = 150nm
pixels = integrate_pixels(psf, camera, emitter)
# Sum of pixels will be ≤ 1000, depending on how much of the PSF is captured by the camera
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
    
    # Get pixel dimensions
    nx = length(camera.pixel_edges_x) - 1
    ny = length(camera.pixel_edges_y) - 1
    
    # Allocate result matrix
    result = Matrix{T}(undef, ny, nx)
    
    # Use in-place function
    integrate_pixels!(
        result,
        psf,
        camera,
        emitter,
        sampling=sampling
    )
    
    return result
end

"""
    integrate_pixels_amplitude(psf::AbstractPSF, camera::AbstractCamera, emitter::AbstractEmitter; sampling::Integer=2)

Integrate PSF complex amplitude over camera pixels.

For each pixel in the camera, numerically integrates the PSF amplitude using the specified
sampling density. Unlike intensity integration, returns unnormalized complex amplitudes
which can be used for coherent calculations.
Automatically handles z-coordinate if both PSF and emitter support it.

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
    
    # Get pixel dimensions
    nx = length(camera.pixel_edges_x) - 1
    ny = length(camera.pixel_edges_y) - 1
    
    # Allocate result matrix
    result = Matrix{T}(undef, ny, nx)
    
    # Use in-place function
    integrate_pixels_amplitude!(
        result,
        psf,
        camera,
        emitter,
        sampling=sampling
    )
    
    return result
end