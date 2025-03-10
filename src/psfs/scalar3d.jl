# src/psfs/scalar3d.jl

"""
    Scalar3DPSF(nₐ::Real, λ::Real, n::Real;
                pupil=nothing,
                pupil_data::Union{Nothing,AbstractMatrix}=nothing,
                coeffs=nothing)

Create a scalar 3D PSF with the specified optical parameters.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `n`: Refractive index of the medium

# Keyword Arguments
- `pupil`: Optional pre-created PupilFunction
- `pupil_data`: Optional complex matrix to initialize the pupil function
- `coeffs`: Optional ZernikeCoefficients for aberrations

# Returns
- `Scalar3DPSF`: 3D PSF with scalar diffraction model

# Notes
- Exactly one of `pupil`, `pupil_data`, or `coeffs` should be provided
- If none are provided, creates an unaberrated pupil
- The constructed PSF stores the provided Zernike coefficients for later use

# Examples
```julia
# Create an unaberrated PSF
psf = Scalar3DPSF(1.4, 0.532, 1.518)

# Create PSF with Zernike aberrations
zc = ZernikeCoefficients(15)
add_astigmatism!(zc, 0.5)  # Add astigmatism
psf = Scalar3DPSF(1.4, 0.532, 1.518; coeffs=zc)
```
"""
function Scalar3DPSF(nₐ::Real, λ::Real, n::Real;
    pupil::Union{Nothing,PupilFunction}=nothing,
    pupil_data::Union{Nothing,AbstractMatrix}=nothing,
    coeffs::Union{Nothing,ZernikeCoefficients}=nothing)
    # Convert inputs to PupilFunction
    pupil_func = if !isnothing(pupil)
        # Direct PupilFunction input
        pupil
    elseif !isnothing(pupil_data)
        # Create PupilFunction from matrix
        PupilFunction(nₐ, λ, n, pupil_data)
    else
        # Default to Zernike representation
        zernike_coeffs = isnothing(coeffs) ? ZernikeCoefficients(1) : coeffs
        PupilFunction(nₐ, λ, n, zernike_coeffs)
    end

    normalize!(pupil_func)
    
    # Store the Zernike coefficients if they were used
    stored_coeffs = !isnothing(coeffs) ? coeffs : nothing
    
    return Scalar3DPSF(nₐ, λ, n, pupil_func, stored_coeffs)
end

"""
    amplitude(psf::Scalar3DPSF{T}, x::Real, y::Real, z::Real) where {T}

Calculate complex amplitude at a 3D position using scalar diffraction theory.

# Arguments
- `psf`: Scalar3DPSF instance
- `x, y, z`: Position in microns relative to PSF center

# Returns
- Complex amplitude at the specified position

# Notes
- Uses Fourier optics to propagate from pupil to image space
- Accounts for defocus and aberrations encoded in the pupil function
"""
function amplitude(psf::Scalar3DPSF{T}, x::Real, y::Real, z::Real) where {T}
    # Use existing PupilFunction amplitude method
    return amplitude(psf.pupil, x, y, z)
end

"""
    (psf::Scalar3DPSF)(x::Real, y::Real, z::Real)

Evaluate PSF intensity at the given 3D position.

# Arguments
- `x, y, z`: Position in microns relative to PSF center

# Returns
- Intensity value at the specified position

# Notes
- Calculated as |amplitude|² of the complex field
- Follows standard operator syntax: `intensity = psf(x, y, z)`
"""
function (psf::Scalar3DPSF)(x::Real, y::Real, z::Real)
    return abs2(amplitude(psf, x, y, z))
end

"""
    integrate_pixels(psf::Scalar3DPSF,
                    camera::AbstractCamera,
                    emitter::AbstractEmitter;
                    sampling::Integer=2)

Integrate PSF intensity over camera pixels for a 3D emitter.

# Arguments
- `psf`: Scalar3DPSF instance
- `camera`: Camera geometry defining pixel layout
- `emitter`: Emitter with position information (must have z-coordinate)
- `sampling`: Subpixel sampling density for integration accuracy (default: 2)

# Returns
- Array of integrated PSF intensities with dimensions [ny, nx]
- Values represent actual photon counts based on emitter.photons

# Notes
- Uses the z-coordinate from the emitter for the focal plane
- Higher sampling values give more accurate results but slower computation
"""
function integrate_pixels(psf::Scalar3DPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2)
    
    # Check if emitter has required z-coordinate
    if !hasfield(typeof(emitter), :z)
        throw(ArgumentError("Scalar3DPSF requires an emitter with a z-coordinate"))
    end
    
    result = _integrate_pixels_generic(psf, camera, emitter,
        (p, x, y) -> p(x, y, emitter.z),
        Float64;
        sampling=sampling)
    
    # Multiply by photon count to preserve physical meaning
    return result .* emitter.photons
end

function Base.show(io::IO, psf::Scalar3DPSF)
    print(io, "Scalar3DPSF(NA=$(psf.nₐ), λ=$(psf.λ)μm, n=$(psf.n))")
    if !isnothing(psf.zernike_coeffs)
        print(io, " with $(length(psf.zernike_coeffs)) Zernike terms")
    end
end