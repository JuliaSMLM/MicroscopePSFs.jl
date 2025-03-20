# src/integration/integration_single.jl

"""
    integrate_pixels!(
        result::AbstractMatrix{T},
        psf::AbstractPSF,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter;
        sampling::Integer=2,
        threaded::Bool=true
    ) where T <: Real

Integrate PSF intensity over camera pixels, storing the result in a pre-allocated matrix.
Automatically uses z-coordinate if both PSF and emitter support it.

# Arguments
- `result`: Pre-allocated array where results will be stored
- `psf`: Point spread function to integrate
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)
  Set to false when using with automatic differentiation frameworks

# Returns
- The `result` array, now filled with integrated intensities scaled by emitter.photons

# Notes
- Array is indexed as [y,x] with [1,1] at top-left pixel
"""
function integrate_pixels!(
    result::AbstractMatrix{T},
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter;
    sampling::Integer=2,
    threaded::Bool=true
) where T <: Real
    
    # Define a function that handles both 2D and 3D cases automatically
    # This is important - we don't explicitly check if we need a 2D or 3D evaluation here
    # Instead, we pass a function that will be determined based on the PSF and emitter traits
    eval_func = if supports_3d(psf) && has_z_coordinate(emitter)
        (p, x, y, z) -> p(x, y, z)
    else
        (p, x, y) -> p(x, y)
    end
    
    # Perform integration
    _integrate_pixels_generic!(
        result,
        psf,
        pixel_edges_x,
        pixel_edges_y,
        emitter,
        eval_func,  # Pass this flexible function
        sampling=sampling,
        threaded=threaded
    )
    
    # Scale by photon count
    result .*= emitter.photons
    
    return result
end

"""
    integrate_pixels!(
        result::AbstractMatrix{T},
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2,
        threaded::Bool=true
    ) where T <: Real

Integrate PSF intensity over camera pixels using camera object, storing the result in a pre-allocated matrix.

# Arguments
- `result`: Pre-allocated array where results will be stored
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- The `result` array, now filled with integrated intensities
"""
function integrate_pixels!(
    result::AbstractMatrix{T},
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2,
    threaded::Bool=true
) where T <: Real
    
    # Delegate to the version that takes pixel edges
    integrate_pixels!(
        result,
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitter;
        sampling=sampling,
        threaded=threaded
    )
end

"""
    integrate_pixels_amplitude!(
        result::AbstractMatrix{Complex{T}},
        psf::AbstractPSF,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter;
        sampling::Integer=2,
        threaded::Bool=true
    ) where T <: Real

Integrate PSF complex amplitude over camera pixels, storing the result in a pre-allocated matrix.
Automatically uses z-coordinate if both PSF and emitter support it.

# Arguments
- `result`: Pre-allocated complex array where results will be stored
- `psf`: Point spread function to integrate
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- The `result` array, now filled with integrated complex amplitudes

# Notes
- Array is indexed as [y,x] with [1,1] at top-left pixel
"""
function integrate_pixels_amplitude!(
    result::AbstractMatrix{Complex{T}},
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter;
    sampling::Integer=2,
    threaded::Bool=true
) where T <: Real
    
    # Use the amplitude function which already has methods for both 2D and 3D
    _integrate_pixels_generic!(
        result,
        psf,
        pixel_edges_x,
        pixel_edges_y,
        emitter,
        amplitude,  # The amplitude function already has appropriate methods
        sampling=sampling,
        threaded=threaded
    )
    
    return result
end

"""
    integrate_pixels_amplitude!(
        result::AbstractMatrix{Complex{T}},
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2,
        threaded::Bool=true
    ) where T <: Real

Integrate PSF complex amplitude over camera pixels using camera object, storing the result in a pre-allocated matrix.

# Arguments
- `result`: Pre-allocated complex array where results will be stored
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- The `result` array, now filled with integrated complex amplitudes
"""
function integrate_pixels_amplitude!(
    result::AbstractMatrix{Complex{T}},
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2,
    threaded::Bool=true
) where T <: Real
    
    # Delegate to the version that takes pixel edges
    integrate_pixels_amplitude!(
        result,
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitter;
        sampling=sampling,
        threaded=threaded
    )
end

"""
    integrate_pixels(
        psf::AbstractPSF, 
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter; 
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF intensity over camera pixels with optional support region optimization.

For each pixel in the specified region, numerically integrates the PSF intensity using the specified 
sampling density. Physical coordinates are relative to camera with (0,0) at top-left corner.
Automatically handles z-coordinate if both PSF and emitter support it.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `pixel_edges_x::AbstractVector`: X-coordinate edges of pixels in microns
- `pixel_edges_y::AbstractVector`: Y-coordinate edges of pixels in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `support`: Region to calculate (default: Inf = full image)
  - If Real: radius in microns around emitter
  - If Tuple: explicit (x_min, x_max, y_min, y_max) in microns
- `sampling::Integer=2`: Number of samples per pixel in each dimension
- `threaded::Bool=true`: Whether to use multi-threading for integration
  Set to false when using with automatic differentiation frameworks

# Returns
- `Matrix{T}`: Integrated intensities where T matches emitter.photons type
- Array is indexed as [y,x] with [1,1] at top-left pixel
- Values represent actual photon counts in each pixel based on emitter's photon value

# Examples
```julia
# Create pixel edges for a 20x20 camera with 100nm pixels
pixel_edges_x = pixel_edges_y = 0:0.1:2.0
emitter = Emitter2D(1.0, 1.0, 1000.0)  # Emitter at (1μm, 1μm) with 1000 photons
psf = GaussianPSF(0.15)  # σ = 150nm

# Calculate over full image
pixels = integrate_pixels(psf, pixel_edges_x, pixel_edges_y, emitter)

# Calculate only within a 0.5μm radius of the emitter
pixels_roi = integrate_pixels(psf, pixel_edges_x, pixel_edges_y, emitter, support=0.5)
```

See also: [`integrate_pixels_amplitude`](@ref), [`AbstractPSF`](@ref)
"""
function integrate_pixels(
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter;
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    T = typeof(emitter.photons)
    
    # Get pixel indices for the support region
    i_range, j_range = get_pixel_indices(pixel_edges_x, pixel_edges_y, emitter, support)
    
    # Create views of the pixel edges for the region
    edges_x_view = view(pixel_edges_x, minimum(i_range):maximum(i_range)+1)
    edges_y_view = view(pixel_edges_y, minimum(j_range):maximum(j_range)+1)
    
    # Allocate result for the full domain
    nx_full = length(pixel_edges_x) - 1
    ny_full = length(pixel_edges_y) - 1
    result_full = zeros(T, ny_full, nx_full)
    
    # Create view of the result for the region
    result_view = view(result_full, j_range, i_range)
    
    # Call in-place function on the views
    integrate_pixels!(
        result_view,
        psf,
        edges_x_view,
        edges_y_view,
        emitter,
        sampling=sampling,
        threaded=threaded
    )
    
    return result_full
end

"""
    integrate_pixels(
        psf::AbstractPSF, 
        camera::AbstractCamera, 
        emitter::AbstractEmitter;
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF intensity over camera pixels with optional support region optimization.

This version takes a camera object instead of explicit pixel edges.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `camera::AbstractCamera`: Camera geometry defining pixel edges in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `support`: Region to calculate (default: Inf = full image)
  - If Real: radius in microns around emitter
  - If Tuple: explicit (x_min, x_max, y_min, y_max) in microns
- `sampling::Integer=2`: Number of samples per pixel in each dimension
- `threaded::Bool=true`: Whether to use multi-threading for integration

# Returns
- Array of pixel values with dimensions [y,x] (matrix indexed [row,col])
- Values are normalized to sum to 1
- Array indices start at [1,1] for top-left pixel

See also: [`integrate_pixels_amplitude`](@ref), [`AbstractPSF`](@ref)
"""
function integrate_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    # Delegate to the version that takes pixel edges
    return integrate_pixels(
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitter;
        support=support,
        sampling=sampling,
        threaded=threaded
    )
end

"""
    integrate_pixels_amplitude(
        psf::AbstractPSF, 
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter; 
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF complex amplitude over camera pixels with optional support region optimization.

For each pixel in the specified region, numerically integrates the PSF amplitude using the specified
sampling density. Unlike intensity integration, returns unnormalized complex amplitudes
which can be used for coherent calculations.
Automatically handles z-coordinate if both PSF and emitter support it.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `pixel_edges_x::AbstractVector`: X-coordinate edges of pixels in microns
- `pixel_edges_y::AbstractVector`: Y-coordinate edges of pixels in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `support`: Region to calculate (default: Inf = full image)
  - If Real: radius in microns around emitter
  - If Tuple: explicit (x_min, x_max, y_min, y_max) in microns
- `sampling::Integer=2`: Number of samples per pixel in each dimension
- `threaded::Bool=true`: Whether to use multi-threading for integration

# Returns
- `Matrix{Complex{T}}`: Integrated complex amplitudes where T matches emitter.photons type
- Array is indexed as [y,x] with [1,1] at top-left pixel
- Values are not normalized to preserve complex amplitude relationships

# Notes
- For coherent calculations, use this function instead of `integrate_pixels`
- Return type is complex to support PSFs with phase information
- To get intensity from amplitude: `abs2.(integrate_pixels_amplitude(...))`

See also: [`integrate_pixels`](@ref), [`AbstractPSF`](@ref)
"""
function integrate_pixels_amplitude(
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitter::AbstractEmitter;
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    T = Complex{typeof(emitter.photons)}
    
    # Get pixel indices for the support region
    i_range, j_range = get_pixel_indices(pixel_edges_x, pixel_edges_y, emitter, support)
    
    # Create views of the pixel edges for the region
    edges_x_view = view(pixel_edges_x, minimum(i_range):maximum(i_range)+1)
    edges_y_view = view(pixel_edges_y, minimum(j_range):maximum(j_range)+1)
    
    # Allocate result for the full domain
    nx_full = length(pixel_edges_x) - 1
    ny_full = length(pixel_edges_y) - 1
    result_full = zeros(T, ny_full, nx_full)
    
    # Create view of the result for the region
    result_view = view(result_full, j_range, i_range)
    
    # Call in-place function on the views
    integrate_pixels_amplitude!(
        result_view,
        psf,
        edges_x_view,
        edges_y_view,
        emitter,
        sampling=sampling,
        threaded=threaded
    )
    
    return result_full
end

"""
    integrate_pixels_amplitude(
        psf::AbstractPSF, 
        camera::AbstractCamera, 
        emitter::AbstractEmitter;
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF complex amplitude over camera pixels with optional support region optimization.

This version takes a camera object instead of explicit pixel edges.

# Arguments
- `psf::AbstractPSF`: Point spread function to integrate
- `camera::AbstractCamera`: Camera geometry defining pixel edges in microns
- `emitter::AbstractEmitter`: Emitter with position in microns relative to camera
- `support`: Region to calculate (default: Inf = full image)
  - If Real: radius in microns around emitter
  - If Tuple: explicit (x_min, x_max, y_min, y_max) in microns
- `sampling::Integer=2`: Number of samples per pixel in each dimension
- `threaded::Bool=true`: Whether to use multi-threading for integration

# Returns
- Matrix of complex amplitudes
- Array indices start at [1,1] for top-left pixel

See also: [`integrate_pixels`](@ref), [`AbstractPSF`](@ref)
"""
function integrate_pixels_amplitude(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    # Delegate to the version that takes pixel edges
    return integrate_pixels_amplitude(
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitter;
        support=support,
        sampling=sampling,
        threaded=threaded
    )
end