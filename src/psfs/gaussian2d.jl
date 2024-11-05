# Type defined in types.jl 

"""
    Gaussian2D{T<:AbstractFloat} <: AbstractPSF

Isotropic 2D Gaussian PSF.

# Fields
- `σ`: Standard deviation in physical units (microns)
"""
struct Gaussian2D{T<:AbstractFloat} <: AbstractPSF
    σ::T

    function Gaussian2D(σ::Real)
        σ > zero(σ) || throw(ArgumentError("σ must be positive"))
        new{float(typeof(σ))}(σ)
    end
end

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

Convert an Airy2D PSF to an approximate Gaussian2D PSF by matching the intensity width.
Uses the empirical relationship that σ ≈ 0.42λ/NA.
"""
function Gaussian2D(psf::Airy2D)
    σ = 0.42 * psf.λ / psf.nₐ
    return Gaussian2D(σ)
end

Base.show(io::IO, psf::Gaussian2D) = print(io, "Gaussian2D(σ=$(psf.σ)μm)")