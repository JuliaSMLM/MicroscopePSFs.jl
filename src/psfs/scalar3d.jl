# Factory function with more flexible input handling
function Scalar3DPSF(nₐ::Real, λ::Real, n::Real; 
                    pupil=nothing, 
                    pupil_data::Union{Nothing, AbstractMatrix}=nothing,
                    coeffs=nothing)
    if !isnothing(pupil)
        # Direct PupilFunction input
        return Scalar3DPupilPSF(nₐ, λ, n, pupil)
    elseif !isnothing(pupil_data)
        # Create PupilFunction from matrix
        pupil_func = PupilFunction(nₐ, λ, n, pupil_data)
        return Scalar3DPupilPSF(nₐ, λ, n, pupil_func)
    else
        # Default to Zernike representation
        coeffs = isnothing(coeffs) ? ZernikeCoefficients(1) : coeffs
        return Scalar3DZernikePSF(nₐ, λ, n, coeffs)
    end
end


function amplitude(psf::Scalar3DZernikePSF{T}, x::Real, y::Real, z::Real) where T
    
    # Generate pupil from Zernike coefficients
    pupil = PupilFunction(psf.nₐ, psf.λ, psf.n, psf.coeffs)

    # Use existing PupilFunction amplitude method
    return amplitude(pupil, x, y, z)
end


function amplitude(psf::Scalar3DPupilPSF{T}, x::Real, y::Real, z::Real) where T
    # Use existing PupilFunction amplitude method
    return amplitude(psf.pupil, x, y, z)
end

# Common intensity method
function (psf::AbstractScalar3DPSF)(x::Real, y::Real, z::Real)
    return abs2(amplitude(psf, x, y, z))
end

# Default pixel integration
function integrate_pixels(psf::AbstractScalar3DPSF, 
                        camera::AbstractCamera, 
                        emitter::AbstractEmitter;
                        sampling::Integer=2)
    return _integrate_pixels_generic(psf, camera, emitter, 
                                   (p,x,y) -> p(x,y,emitter.z), 
                                   Float64;
                                   sampling=sampling)
end