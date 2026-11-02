# src/psfs/airy2d.jl 

"""
    (psf::AiryPSF)(x, y)

Airy pattern intensity for circular aperture under paraxial approximation.
"""
function (psf::AiryPSF)(x::Real, y::Real)
    a = amplitude(psf, x, y)
    return abs2(a)  # Could also write this calculation out explicitly if we want
end

function amplitude(psf::AiryPSF, x::Real, y::Real)
    r² = x^2 + y^2
    w² = r² * psf.ν^2
    # Series of 2J₁(w)/w; no square root at r = 0, so derivatives survive.
    if w² < 1.0e-10
        return psf.ν / sqrt(4π) * (1 - w² / 8)
    end
    w = sqrt(r²) * psf.ν
    return psf.ν / sqrt(4π) * (2 * besselj1(w) / w)
end



"""
    AiryPSF(psf::GaussianPSF, λ)

Convert a GaussianPSF to an approximate AiryPSF.

Uses the inverse of the Gaussian conversion relationship.

"""
function AiryPSF(psf::GaussianPSF; λ::Real=0.532)
    # Invert the relationship σ = 0.42λ/NA
    nₐ = 0.22 * λ / psf.σ
    return AiryPSF(nₐ, λ)
end



# Add pretty printing
function Base.show(io::IO, psf::AiryPSF)
    print(io, "AiryPSF(NA=$(psf.nₐ), λ=$(psf.λ)μm)")
end

