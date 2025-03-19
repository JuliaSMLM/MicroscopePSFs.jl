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
    GaussianPSF{T<:AbstractFloat} <: Abstract2DPSF{T}

Isotropic 2D Gaussian PSF.

# Fields
- `σ`: Standard deviation in physical units (typically microns)
"""
struct GaussianPSF{T<:AbstractFloat} <: Abstract2DPSF{T}
    σ::T
end

"""
    AiryPSF{T<:AbstractFloat} <: Abstract2DPSF{T}

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
struct AiryPSF{T<:AbstractFloat} <: Abstract2DPSF{T}
    nₐ::T
    λ::T
    ν::T

    function AiryPSF(nₐ::Real, λ::Real)
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
    ScalarPSF{T} <: Abstract3DPSF{T}

Scalar 3D PSF using explicit pupil function representation.

# Fields
- `nₐ::T`: Numerical aperture
- `λ::T`: Wavelength in microns
- `n::T`: Refractive index
- `pupil::PupilFunction{T}`: Complex pupil function
- `zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}`: Zernike coefficients used to create this PSF (if applicable)

# Notes
Can be initialized with either a PupilFunction or ZernikeCoefficients
using the ScalarPSF factory function.
"""
struct ScalarPSF{T} <: Abstract3DPSF{T}
    nₐ::T
    λ::T
    n::T
    pupil::PupilFunction{T}
    zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}

    function ScalarPSF(nₐ::Real, λ::Real, n::Real, pupil::PupilFunction, zernike_coeffs::Union{Nothing, ZernikeCoefficients}=nothing)
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), real(eltype(pupil.field)))
        nₐ > 0 || throw(ArgumentError("NA must be positive"))
        λ > 0 || throw(ArgumentError("wavelength must be positive"))
        n > 0 || throw(ArgumentError("refractive index must be positive"))
        new{T}(T(nₐ), T(λ), T(n), pupil, zernike_coeffs)
    end
end

"""
    VectorPSF{T<:AbstractFloat} <: Abstract3DPSF{T}

Vector PSF model using explicit pupil function representation.
Allows direct manipulation of pupil function for custom aberrations.

# Fields
- `nₐ::T`: Numerical aperture
- `λ::T`: Wavelength in microns
- `n_medium::T`: Sample medium refractive index
- `n_coverslip::T`: Cover slip refractive index
- `n_immersion::T`: Immersion medium refractive index
- `dipole::DipoleVector{T}`: Dipole orientation
- `z_stage::T`: Distance the sample stage was moved away from the nominal focal plane at the coverslip (μm)
- `vector_pupils::Vector{VectorPupilFunction{T}}`: Vector of pupil functions containing field components.
  For a single dipole, this contains one pupil. For a rotating dipole, it contains three pupils
  for x, y, and z orientations.
- `base_pupil::Union{Nothing, PupilFunction{T}}`: Base pupil function representing system aberrations
- `zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}`: Zernike coefficients used to create this PSF (if applicable)
"""
struct VectorPSF{T<:AbstractFloat} <: Abstract3DPSF{T}
    nₐ::T
    λ::T
    n_medium::T
    n_coverslip::T
    n_immersion::T
    dipole::DipoleVector{T}
    z_stage::T
    vector_pupils::Vector{VectorPupilFunction{T}}
    base_pupil::Union{Nothing, PupilFunction{T}}
    zernike_coeffs::Union{Nothing, ZernikeCoefficients{T}}
end

"""
    SplinePSF{T<:AbstractFloat, IT<:AbstractInterpolation} <: AbstractPSF

A point spread function (PSF) represented as a B-spline interpolation.

# Fields
- `spline`: The B-spline interpolation object 
- `x_range`: Range of x-coordinates used for uniform grid interpolation
- `y_range`: Range of y-coordinates used for uniform grid interpolation  
- `z_range`: Range of z-coordinates for 3D PSFs, or `nothing` for 2D PSFs
- `original_grid`: Original grid data used to create the interpolation
- `interp_order`: Interpolation order used (0=constant, 1=linear, 3=cubic)

# Notes
- Coordinates and ranges are in physical units (typically microns)
- PSF values are preserved from the original PSF that was sampled
- Full implementation is in spline_psf.jl
"""
struct SplinePSF{T<:AbstractFloat, IT<:AbstractInterpolation} <: AbstractPSF
    spline::IT
    x_range::StepRangeLen{T}
    y_range::StepRangeLen{T}
    z_range::Union{StepRangeLen{T}, Nothing}
    original_grid::Array{T}
    interp_order::Int
end

