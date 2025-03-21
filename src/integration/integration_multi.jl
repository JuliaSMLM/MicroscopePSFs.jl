# src/integration/integration_multi.jl

"""
    integrate_pixels(
        psf::AbstractPSF,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitters::Vector{<:AbstractEmitter};
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF intensity over camera pixels for multiple emitters with optional support region optimization.

# Arguments
- `psf`: Point spread function to integrate
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `emitters`: Vector of emitters with position information
- `support`: Region to calculate for each emitter (default: Inf = full image)
  - If Real: radius in microns around each emitter
  - If Tuple: explicit (x_min, x_max, y_min, y_max) in microns (fixed region for all emitters)
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)
  Set to false when using with automatic differentiation frameworks

# Returns
- Array of integrated PSF intensities with dimensions matching the camera
- Values represent actual photon counts from all emitters

# Notes
- Results are the sum of individual emitter contributions (incoherent addition)
- For coherent addition, use `integrate_pixels_amplitude` and sum complex amplitudes
"""
function integrate_pixels(
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitters::Vector{<:AbstractEmitter};
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    # Return zeros if no emitters
    isempty(emitters) && return zeros(Float64, length(pixel_edges_y)-1, length(pixel_edges_x)-1)
    
    # Determine result type from first emitter (promotion will happen automatically)
    T = typeof(emitters[1].photons)
    
    # Allocate result array for the full domain
    nx = length(pixel_edges_x) - 1
    ny = length(pixel_edges_y) - 1
    result = zeros(T, ny, nx)
    
    # Process each emitter by calling the single emitter version and accumulating results
    for emitter in emitters
        # Get the contribution from this emitter
        emitter_result = integrate_pixels(
            psf,
            pixel_edges_x,
            pixel_edges_y,
            emitter;
            support=support,
            sampling=sampling,
            threaded=threaded
        )
        
        # Add to the result
        result .+= emitter_result
    end
    
    return result
end

"""
    integrate_pixels(
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitters::Vector{<:AbstractEmitter};
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF intensity over camera pixels for multiple emitters with optional support region optimization.
This version takes a camera object instead of explicit pixel edges.

# Arguments
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitters`: Vector of emitters with position information
- `support`: Region to calculate for each emitter (default: Inf = full image)
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- Array of integrated PSF intensities with dimensions matching the camera
- Values represent actual photon counts from all emitters
"""
function integrate_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitters::Vector{<:AbstractEmitter};
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    # Delegate to the version that takes pixel edges
    return integrate_pixels(
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitters;
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
        emitters::Vector{<:AbstractEmitter};
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF complex amplitude over camera pixels for multiple emitters with optional support region optimization.

# Arguments
- `psf`: Point spread function to integrate
- `pixel_edges_x`, `pixel_edges_y`: Arrays defining pixel edge coordinates
- `emitters`: Vector of emitters with position information
- `support`: Region to calculate for each emitter (default: Inf = full image)
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- Array of integrated complex amplitudes
- Values represent coherently summed field contributions

# Notes
- Results are the coherent sum of individual emitter field contributions
- For incoherent addition, use `abs2.(integrate_pixels_amplitude(...))` or use `integrate_pixels`
"""
function integrate_pixels_amplitude(
    psf::AbstractPSF,
    pixel_edges_x::AbstractVector,
    pixel_edges_y::AbstractVector,
    emitters::Vector{<:AbstractEmitter};
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    # Return zeros if no emitters
    isempty(emitters) && return zeros(Complex{Float64}, length(pixel_edges_y)-1, length(pixel_edges_x)-1)
    
    # Determine result type from first emitter
    T = Complex{typeof(emitters[1].photons)}
    
    # Allocate result array for the full domain
    nx = length(pixel_edges_x) - 1
    ny = length(pixel_edges_y) - 1
    result = zeros(T, ny, nx)
    
    # Process each emitter by calling the single emitter version and accumulating results
    for emitter in emitters
        # Get the contribution from this emitter
        emitter_result = integrate_pixels_amplitude(
            psf,
            pixel_edges_x,
            pixel_edges_y,
            emitter;
            support=support,
            sampling=sampling,
            threaded=threaded
        )
        
        # Add to the result
        result .+= emitter_result
    end
    
    return result
end

"""
    integrate_pixels_amplitude(
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitters::Vector{<:AbstractEmitter};
        support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF complex amplitude over camera pixels for multiple emitters with optional support region optimization.
This version takes a camera object instead of explicit pixel edges.

# Arguments
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitters`: Vector of emitters with position information
- `support`: Region to calculate for each emitter (default: Inf = full image)
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- Array of integrated complex amplitudes
- Values represent coherently summed field contributions
"""
function integrate_pixels_amplitude(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitters::Vector{<:AbstractEmitter};
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}} = Inf,
    sampling::Integer=2,
    threaded::Bool=true
)
    # Delegate to the version that takes pixel edges
    return integrate_pixels_amplitude(
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitters;
        support=support,
        sampling=sampling,
        threaded=threaded
    )
end