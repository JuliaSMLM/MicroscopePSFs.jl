# src/integration/integration_multi.jl

"""
    integrate_pixels(
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitters::Vector{<:AbstractEmitter};
        sampling::Integer=2
    )

Integrate PSF intensity over camera pixels for multiple emitters.

# Arguments
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitters`: Vector of emitters with position information
- `sampling`: Subpixel sampling density (default: 2)

# Returns
- Array of integrated PSF intensities with dimensions matching the camera
- Values represent actual photon counts from all emitters

# Notes
- Results are the sum of individual emitter contributions (incoherent addition)
- For coherent addition, use `integrate_pixels_amplitude` and sum complex amplitudes
"""
function integrate_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitters::Vector{<:AbstractEmitter};
    sampling::Integer=2
)
    isempty(emitters) && return zeros(Float64, length(camera.pixel_edges_y)-1, length(camera.pixel_edges_x)-1)
    
    # Determine result type from first emitter (promotion will happen automatically)
    T = typeof(emitters[1].photons)
    
    # Allocate result array
    ny = length(camera.pixel_edges_y) - 1
    nx = length(camera.pixel_edges_x) - 1
    result = zeros(T, ny, nx)
    
    # Temporary buffer for individual emitter contribution
    buffer = similar(result)
    
    # Process each emitter
    for emitter in emitters
        fill!(buffer, zero(T))
        integrate_pixels!(buffer, psf, camera, emitter; sampling=sampling)
        result .+= buffer
    end
    
    return result
end

"""
    integrate_pixels_amplitude(
        psf::AbstractPSF,
        camera::AbstractCamera,
        emitters::Vector{<:AbstractEmitter};
        sampling::Integer=2
    )

Integrate PSF complex amplitude over camera pixels for multiple emitters.

# Arguments
- `psf`: Point spread function to integrate
- `camera`: Camera geometry defining pixel edges
- `emitters`: Vector of emitters with position information
- `sampling`: Subpixel sampling density (default: 2)

# Returns
- Array of integrated complex amplitudes
- Values represent coherently summed field contributions

# Notes
- Results are the coherent sum of individual emitter field contributions
- For incoherent addition, use `abs2.(integrate_pixels_amplitude(...))` or use `integrate_pixels`
"""
function integrate_pixels_amplitude(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitters::Vector{<:AbstractEmitter};
    sampling::Integer=2
)
    isempty(emitters) && return zeros(Complex{Float64}, length(camera.pixel_edges_y)-1, length(camera.pixel_edges_x)-1)
    
    # Determine result type from first emitter
    T = Complex{typeof(emitters[1].photons)}
    
    # Allocate result array
    ny = length(camera.pixel_edges_y) - 1
    nx = length(camera.pixel_edges_x) - 1
    result = zeros(T, ny, nx)
    
    # Temporary buffer for individual emitter contribution
    buffer = similar(result)
    
    # Process each emitter
    for emitter in emitters
        fill!(buffer, zero(T))
        integrate_pixels_amplitude!(buffer, psf, camera, emitter; sampling=sampling)
        result .+= buffer
    end
    
    return result
end