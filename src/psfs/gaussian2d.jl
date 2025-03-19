# src/psfs/gaussian2d.jl

# Type defined in types.jl 

"""
    (psf::GaussianPSF)(x, y)

Gaussian PSF intensity. Radially symmetric with standard deviation σ.
"""
function (psf::GaussianPSF)(x::Real, y::Real)
    r² = x^2 + y^2
    return 1/(2π*psf.σ^2) * exp(-r²/(2psf.σ^2))
end


function amplitude(psf::GaussianPSF, x::Real, y::Real)
    return sqrt(psf(x, y))
end




"""
    GaussianPSF(psf::AiryPSF)

Convert an AiryPSF to an approximate GaussianPSF by matching the FWHM.

Uses the relationship that σ ≈ 0.22 * λ/NA
"""
function GaussianPSF(psf::AiryPSF)
    σ = 0.22 * psf.λ / psf.nₐ
    return GaussianPSF(σ)
end

Base.show(io::IO, psf::GaussianPSF) = print(io, "GaussianPSF(σ=$(psf.σ)μm)")