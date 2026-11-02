# src/pupil.jl

"""
    PupilFunction{T}

Represents a pupil function with physical parameters.

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in μm
- `n`: Refractive index
- `field`: Square complex field array with at least two samples per axis
"""
struct PupilFunction{T<:AbstractFloat}
    nₐ::T
    λ::T
    n::T
    field::Array{Complex{T},2}

    function PupilFunction(nₐ::Real, λ::Real, n::Real, field::Array{<:Complex})
        nₐ > 0 || throw(ArgumentError("Numerical aperture must be positive"))
        λ > 0 || throw(ArgumentError("Wavelength must be positive"))
        n > 0 || throw(ArgumentError("Refractive index must be positive"))
        ndims(field) == 2 && size(field, 1) == size(field, 2) &&
            size(field, 1) >= 2 || throw(ArgumentError(
                "Pupil field must be square with at least two samples per axis; " *
                "supply a square field"))
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), real(eltype(field)))
        new{T}(T(nₐ), T(λ), T(n), convert(Array{Complex{T}}, field))
    end
end

# Computed properties
"""
    kmax(p::PupilFunction)

Get maximum spatial frequency in μm⁻¹"""
kmax(p::PupilFunction) = p.nₐ / p.λ

"""
    k₀(p::PupilFunction)

Get central wavevector magnitude in μm⁻¹"""
k₀(p::PupilFunction) = p.n / p.λ

"""
    kpixelsize(p::PupilFunction)

Get pupil plane sampling in μm⁻¹"""
kpixelsize(p::PupilFunction) = 2kmax(p) / (size(p.field, 1) - 1)

"""
    copy(p::PupilFunction)

Return a deep copy of `p` with an independent `field` array.
"""
Base.copy(p::PupilFunction) = PupilFunction(p.nₐ, p.λ, p.n, copy(p.field))

"""
    PupilFunction(nₐ::Real, λ::Real, n::Real, 
                  zc::ZernikeCoefficients;
                  grid_size::Int=64)

Create a PupilFunction from ZernikeCoefficients type.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `n`: Refractive index
- `zc`: ZernikeCoefficients containing amplitude and phase coefficients

# Keyword Arguments
- `grid_size`: Size of square sampling grid (default: 64). The image intensity repeats
  in x and y with period `L = (grid_size - 1) * λ / (2nₐ)` in microns. Keep render
  `support` below `L/2` or the frame narrower than `L`. A larger grid widens `L`.
  For `grid_size=128`, NA 1.4 and λ 0.6 microns, `L ≈ 27.21` microns.

# Returns
- `PupilFunction` with complex field computed from Zernike polynomials

Amplitude is the sum of the magnitude coefficients times their Zernike modes.
The default unit piston supplies the unaberrated amplitude; no offset is added.
"""
function PupilFunction(nₐ::Real, λ::Real, n::Real,
    zc::ZernikeCoefficients;
    grid_size::Int=64)
    T = promote_type(typeof(nₐ), typeof(λ), typeof(n), eltype(zc.phase))
    field = Zernike._evaluate_pupil(zc, grid_size, T)

    return PupilFunction(nₐ, λ, n, field)
end


# Energy uses exactly the same grid and strict aperture as amplitude evaluation.
function _pupil_energy(p::PupilFunction{T}) where {T}
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    center = (sz + 1) / 2
    energy = zero(T)
    for i in 1:sz, j in 1:sz
        kx = kpix * (i - center)
        ky = kpix * (j - center)
        if kx^2 + ky^2 < kmax(p)^2
            energy += abs2(p.field[j, i])
        end
    end
    return energy * kpix^2
end

# Normalization
"""
    normalize!(p::PupilFunction)

Normalize to unit energy inside the strict circular aperture used by `amplitude`.
The field must have finite positive energy there; zero fields may be stored but
cannot be normalized. Samples outside the aperture do not contribute.
"""
function normalize!(p::PupilFunction)
    # Energy normalization using Parseval's theorem
    total_energy = _pupil_energy(p)
    isfinite(total_energy) && total_energy > 0 || throw(ArgumentError(
        "Pupil energy must be finite and positive inside the aperture; " *
        "supply a finite nonzero field"))
    p.field ./= sqrt(total_energy)
    return p
end

# Display and visualization
function Base.show(io::IO, p::PupilFunction)
    sz = size(p.field, 1)
    print(io, "PupilFunction(NA=$(p.nₐ), λ=$(p.λ)μm, n=$(p.n), $(sz)x$(sz))")
end

# Utility functions for pupil manipulation
"""
    apply_defocus!(p::PupilFunction, z::Real)

Apply defocus phase to pupil function for propagation distance z.
"""
function apply_defocus!(p::PupilFunction, z::Real)
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    k0_center = (sz + 1) / 2

    for i in 1:sz, j in 1:sz
        kx = kpix * (i - k0_center)
        ky = kpix * (j - k0_center)
        kr2 = kx^2 + ky^2

        if kr2 < kmax(p)^2
            kz = sqrt(complex(k₀(p)^2 - kr2))
            p.field[j, i] *= exp(2π * im * z * kz)
        end
    end
    return p
end

"""
    apply_aperture!(p::PupilFunction, radius::Real=1.0)

Apply circular aperture to pupil function. Radius is relative to NA.
"""
function apply_aperture!(p::PupilFunction, radius::Real=1.0)
    sz = size(p.field, 1)
    xs = ys = range(-1, 1, length=sz)

    for i in 1:sz, j in 1:sz
        x, y = xs[i], ys[j]
        if sqrt(x^2 + y^2) > radius
            p.field[j, i] = 0
        end
    end
    return p
end

"""
    amplitude(p::PupilFunction, x::Real, y::Real, z::Real)

Calculate complex amplitude from pupil function integration.

# Arguments
- `p::PupilFunction`: Pupil function
- `x::Real`: X position in μm
- `y::Real`: Y position in μm
- `z::Real`: Z position in μm"""
function amplitude(p::PupilFunction{T}, x::Real, y::Real, z::Real) where {T}
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    # k0_center = (sz + 1) ÷ 2
    k0_center = (sz + 1) / 2

    result = zero(Complex{T})
    kmax² = kmax(p)^2
    k₀² = k₀(p)^2

    @inbounds for i in 1:sz, j in 1:sz
        kx = kpix * (i - k0_center)
        ky = kpix * (j - k0_center)
        kr2 = kx^2 + ky^2

        if kr2 < kmax²
            kz = sqrt(complex(k₀² - kr2))
            phase = 2π * (x * kx + y * ky + z * kz)
            result += p.field[j, i] * exp(im * phase)
        end
    end

    return result * kpix^2
end
