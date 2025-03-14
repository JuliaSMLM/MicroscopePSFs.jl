# src/psfs/gaussian2d.jl

# Type defined in types.jl 

"""
    (psf::Gaussian2D)(x, y)

Gaussian PSF intensity. Radially symmetric with standard deviation σ.
"""
function (psf::Gaussian2D)(x::Real, y::Real)
    r² = x^2 + y^2
    return 1/(2π*psf.σ^2) * exp(-r²/(2psf.σ^2))
end


function amplitude(psf::Gaussian2D, x::Real, y::Real)
    return sqrt(psf(x, y))
end




"""
    Gaussian2D(psf::Airy2D)

Convert an Airy2D PSF to an approximate Gaussian2D PSF by matching the FWHM.

Uses the relationship that σ ≈ 0.22 * λ/NA
"""
function Gaussian2D(psf::Airy2D)
    σ = 0.22 * psf.λ / psf.nₐ
    return Gaussian2D(σ)
end

Base.show(io::IO, psf::Gaussian2D) = print(io, "Gaussian2D(σ=$(psf.σ)μm)")