# src/psfs/types.jl

"""
    AbstractPSF

Abstract base type for all Point Spread Functions
"""
abstract type AbstractPSF end

# 2D PSF Types
"""
    Abstract2DPSF{T<:AbstractFloat}

Abstract type for all 2D point spread functions.
Parameterized by numeric type T.
"""
abstract type Abstract2DPSF{T<:AbstractFloat} <: AbstractPSF end

"""
    Gaussian2D{T<:AbstractFloat} <: Abstract2DPSF{T}

Isotropic 2D Gaussian PSF.

# Fields
- `σ`: Standard deviation in physical units (typically microns)
"""
struct Gaussian2D{T<:AbstractFloat} <: Abstract2DPSF{T}
    σ::T
end

"""
    Airy2D{T<:AbstractFloat} <: Abstract2DPSF{T}

2D Airy pattern PSF using paraxial, scalar model.

Field amplitude for a circular aperture is given by:
    A(r) = ν/√(4π) * (2J₁(νr)/(νr))
where:
    ν = 2π*nₐ/λ
    J₁ is the Bessel function of first kind, order 1
    r is the radial distance from optical axis

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `ν`: Optical parameter = 2π*nₐ/λ
"""
struct Airy2D{T<:AbstractFloat} <: Abstract2DPSF{T}
    nₐ::T
    λ::T
    ν::T

    function Airy2D(nₐ::Real, λ::Real)
        T = promote_type(typeof(nₐ), typeof(λ))
        nₐ > 0 || throw(ArgumentError("Numerical aperture must be positive"))
        λ > 0 || throw(ArgumentError("Wavelength must be positive"))
        ν = 2T(π) * nₐ / λ
        new{T}(T(nₐ), T(λ), ν)
    end
end

# 3D PSF Types
"""
    Abstract3DPSF{T<:AbstractFloat}

Abstract type for all 3D point spread functions.
"""
abstract type Abstract3DPSF{T<:AbstractFloat} <: AbstractPSF end

"""
    Scalar3DPSF{T} <: Abstract3DPSF{T}

Scalar 3D PSF using explicit pupil function representation.

# Fields
- `nₐ::T`: Numerical aperture
- `λ::T`: Wavelength in microns
- `n::T`: Refractive index
- `pupil::PupilFunction{T}`: Complex pupil function
- `zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}`: Zernike coefficients used to create this PSF (if applicable)

# Notes
Can be initialized with either a PupilFunction or ZernikeCoefficients
using the Scalar3DPSF factory function.
"""
struct Scalar3DPSF{T} <: Abstract3DPSF{T}
    nₐ::T
    λ::T
    n::T
    pupil::PupilFunction{T}
    zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}

    function Scalar3DPSF(nₐ::Real, λ::Real, n::Real, pupil::PupilFunction, zernike_coeffs::Union{Nothing, ZernikeCoefficients}=nothing)
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), real(eltype(pupil.field)))
        nₐ > 0 || throw(ArgumentError("NA must be positive"))
        λ > 0 || throw(ArgumentError("wavelength must be positive"))
        n > 0 || throw(ArgumentError("refractive index must be positive"))
        new{T}(T(nₐ), T(λ), T(n), pupil, zernike_coeffs)
    end
end

"""
    Vector3DPSF{T<:AbstractFloat} <: Abstract3DPSF{T}

Vector PSF model using explicit pupil function representation.
Allows direct manipulation of pupil function for custom aberrations.

# Fields
- `nₐ::T`: Numerical aperture
- `λ::T`: Wavelength in microns
- `n_medium::T`: Sample medium refractive index
- `n_coverslip::T`: Cover slip refractive index
- `n_immersion::T`: Immersion medium refractive index
- `dipole::DipoleVector{T}`: Dipole orientation
- `focal_z::T`: Focal plane position (z) in microns relative to nominal focus
- `vector_pupils::VectorPupilFunction{T}`: Pre-calculated pupil functions containing vector field components (Ex,Ey),
  dipole orientation effects, base aberrations, apodization, and all position-independent factors
- `base_pupil::Union{Nothing, PupilFunction{T}}`: Base pupil function representing system aberrations
- `zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}`: Zernike coefficients used to create this PSF (if applicable)
"""
struct Vector3DPSF{T<:AbstractFloat} <: Abstract3DPSF{T}
    nₐ::T
    λ::T
    n_medium::T
    n_coverslip::T
    n_immersion::T
    dipole::DipoleVector{T}
    focal_z::T
    vector_pupils::VectorPupilFunction{T}
    base_pupil::Union{Nothing, PupilFunction{T}}
    zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}
end