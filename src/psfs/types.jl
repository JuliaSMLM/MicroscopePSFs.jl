# Define types here to avoid circular dependencies on conversion functions


"""
    AbstractPSF

Abstract base type for all Point Spread Functions
"""
abstract type AbstractPSF end

"""
    Gaussian2D{T<:AbstractFloat} <: AbstractPSF

Isotropic 2D Gaussian PSF.

# Fields
- `σ`: Standard deviation in physical units (typically microns)
"""
struct Gaussian2D{T<:AbstractFloat} <: AbstractPSF
    σ::T
end

"""
    Airy2D{T<:AbstractFloat} <: AbstractPSF

2D Airy pattern PSF using paraxial, scalar model.

Field amplitude for a circular aperture is given by:
    A(r) = ν/√(4π) * (2J₁(νr)/(νr))
where:
    ν = 2π*nₐ/λ
    J₁ is the Bessel function of first kind, order 1
    r is the radial distance from optical axis

# Fields
- `nₐ`: Numerical Aperture
- `λ`: Wavelength in microns
- `ν`: Optical parameter = 2π*nₐ/λ

Physical coordinates are in microns relative to PSF center.
"""
struct Airy2D{T<:AbstractFloat} <: AbstractPSF
    nₐ::T
    λ::T 
    ν::T

    function Airy2D(nₐ::Real, λ::Real)
        T = promote_type(typeof(nₐ), typeof(λ))
        nₐ > 0 || throw(ArgumentError("Numerical aperture must be positive"))
        λ > 0 || throw(ArgumentError("Wavelength must be positive"))
        ν = 2T(π)*nₐ/λ
        new{T}(T(nₐ), T(λ), ν)    
    end
end
