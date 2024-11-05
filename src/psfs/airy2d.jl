"""
    (psf::Airy2D)(x, y)

Airy pattern intensity for circular aperture under paraxial approximation.
"""
function (psf::Airy2D)(x::Real, y::Real)
    a = amplitude(psf, x, y)
    return abs2(a)  # Could also write this calculation out explicitly if we want
end

function amplitude(psf::Airy2D, x::Real, y::Real)
    r = sqrt(x^2 + y^2)
    w = r * psf.ν
    if w < 1e-5  # Handle center point using series expansion
        return psf.ν/sqrt(4π)  # First term of 2J₁(w)/w series is 1
    else
        return psf.ν/sqrt(4π) * (2*besselj1(w)/w)
    end
end



"""
    Airy2D(psf::Gaussian2D; λ=0.532)

Convert a Gaussian2D PSF to an approximate Airy2D PSF.
Uses the inverse of the Gaussian conversion relationship.
Optional wavelength parameter defaults to 532 nm (green light).

Note: Since Gaussian2D doesn't contain wavelength information,
you may want to specify the wavelength explicitly for more accurate conversion.
"""
function Airy2D(psf::Gaussian2D; λ::Real=0.532)
    # Invert the relationship σ = 0.42λ/NA
    nₐ = 0.42 * λ / psf.σ
    return Airy2D(nₐ, λ)
end



# Add pretty printing
function Base.show(io::IO, psf::Airy2D)
    print(io, "Airy2D(NA=$(psf.nₐ), λ=$(psf.λ)μm)")
end

