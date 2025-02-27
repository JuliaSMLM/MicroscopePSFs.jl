# src/psfs/scalar3d.jl

# Factory function with more flexible input handling
function Scalar3DPSF(nₐ::Real, λ::Real, n::Real;
    pupil=nothing,
    pupil_data::Union{Nothing,AbstractMatrix}=nothing,
    coeffs=nothing)
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
    return Scalar3DPSF(nₐ, λ, n, pupil_func)
end

function amplitude(psf::Scalar3DPSF{T}, x::Real, y::Real, z::Real) where {T}
    # Use existing PupilFunction amplitude method
    return amplitude(psf.pupil, x, y, z)
end

# Intensity method
function (psf::Scalar3DPSF)(x::Real, y::Real, z::Real)
    return abs2(amplitude(psf, x, y, z))
end

# Default pixel integration
function integrate_pixels(psf::Scalar3DPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2)
    
    result = _integrate_pixels_generic(psf, camera, emitter,
        (p, x, y) -> p(x, y, emitter.z),
        Float64;
        sampling=sampling)
    
    # Multiply by photon count to preserve physical meaning
    return result .* emitter.photons
end

# Display method
function Base.show(io::IO, psf::Scalar3DPSF)
    print(io, "Scalar3DPSF(NA=$(psf.nₐ), λ=$(psf.λ)μm, n=$(psf.n))")
end