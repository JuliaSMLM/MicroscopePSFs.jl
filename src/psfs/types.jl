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
        ν = 2T(π)*nₐ/λ
        new{T}(T(nₐ), T(λ), ν)    
    end
end

# 3D PSF Types
"""
    Abstract3DPSF{T<:AbstractFloat}

Abstract type for all 3D point spread functions.
Parameterized by numeric type T.
"""
abstract type Abstract3DPSF{T<:AbstractFloat} <: AbstractPSF end

"""
    AbstractScalar3DPSF{T<:AbstractFloat}

Abstract type for scalar 3D point spread functions.
Implements scalar diffraction theory.
"""
abstract type AbstractScalar3DPSF{T} <: Abstract3DPSF{T} end

"""
    AbstractVector3DPSF{T<:AbstractFloat}

Abstract type for vectorial 3D point spread functions.
Implements full electromagnetic theory.
"""
abstract type AbstractVector3DPSF{T} <: Abstract3DPSF{T} end

"""
    Scalar3DZernikePSF{T<:AbstractFloat}

Scalar 3D PSF using Zernike polynomial representation of pupil function.

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `n`: Refractive index
- `coeffs`: Zernike coefficients for pupil function
"""
struct Scalar3DZernikePSF{T} <: AbstractScalar3DPSF{T}
    nₐ::T
    λ::T
    n::T
    coeffs::ZernikeCoefficients{T}

    function Scalar3DZernikePSF(nₐ::Real, λ::Real, n::Real, coeffs::ZernikeCoefficients)
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), eltype(coeffs.mag))
        nₐ > 0 || throw(ArgumentError("NA must be positive"))
        λ > 0 || throw(ArgumentError("wavelength must be positive"))
        n > 0 || throw(ArgumentError("refractive index must be positive"))
        new{T}(T(nₐ), T(λ), T(n), coeffs)
    end
end

"""
    Scalar3DPupilPSF{T<:AbstractFloat}

Scalar 3D PSF using direct pupil function representation.

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `n`: Refractive index
- `pupil`: Complex pupil function
"""
struct Scalar3DPupilPSF{T} <: AbstractScalar3DPSF{T}
    nₐ::T
    λ::T
    n::T
    pupil::PupilFunction{T}

    function Scalar3DPupilPSF(nₐ::Real, λ::Real, n::Real, pupil::PupilFunction)
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), real(eltype(pupil.field)))
        nₐ > 0 || throw(ArgumentError("NA must be positive"))
        λ > 0 || throw(ArgumentError("wavelength must be positive"))
        n > 0 || throw(ArgumentError("refractive index must be positive"))
        new{T}(T(nₐ), T(λ), T(n), pupil)
    end
end