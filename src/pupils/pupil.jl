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

# Half width of the sampled pupil's alias-free image field, in emitter coordinates.
_alias_free_halfwidth(p) = inv(2kpixelsize(p))

# AD extensions unwrap coordinates only for window membership, keeping their precision.
_primal(x::Real) = x

@inline function _outside_field_window(x::Real, halfwidth::Real)
    p = _primal(x)
    return p <= -halfwidth || p >= halfwidth
end

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
- `grid_size`: Size of square sampling grid (default: 64). The field is zero at
  and beyond `L/2` in either lateral coordinate relative to the emitter, where
  `L = (grid_size - 1) * λ / (2nₐ)` in microns. A larger `grid_size` widens the field.
  For grid 128, NA 1.4 and λ 0.6 microns, `L/2 ≈ 13.61` microns.

# Returns
- `PupilFunction` with complex field computed from Zernike polynomials

Magnitude coefficients are absolute. The amplitude is `Σⱼ zc.mag[j] Zⱼ` with no offset
added, so unit amplitude needs `zc.mag[1] = 1`, which `ZernikeCoefficients(n)` sets.
A caller holding offsets from unit amplitude must add 1 to `zc.mag[1]`: `1 + 0.06 Z₅` is
`mag[1] = 1, mag[5] = 0.06`, while `mag[5] = 0.06` alone gives `0.06 Z₅` and no error.
With every magnitude zero the field is zero, so a phase-only pupil also needs
`mag[1] = 1`; `ScalarPSF(...; zernike_coeffs)` and `VectorPSF(...; base_zernike)` refuse
that field.
"""
function PupilFunction(nₐ::Real, λ::Real, n::Real,
    zc::ZernikeCoefficients;
    grid_size::Int=64)
    T = promote_type(typeof(nₐ), typeof(λ), typeof(n), eltype(zc.phase))
    field = Zernike._evaluate_pupil(zc, grid_size, T)

    return PupilFunction(nₐ, λ, n, field)
end


# Aperture membership uses the same grid arithmetic for every evaluation precision.
@inline function _pupil_sample(p::PupilFunction, i, j)
    center = (size(p.field, 1) + 1) / 2
    kpix = kpixelsize(p)
    kx = kpix * (i - center)
    ky = kpix * (j - center)
    kr2 = kx^2 + ky^2
    return kx, ky, kr2, kr2 < kmax(p)^2
end

function _pupil_energy(p::PupilFunction{T}) where {T}
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    energy = zero(T)
    for i in 1:sz, j in 1:sz
        _, _, _, inside = _pupil_sample(p, i, j)
        if inside
            energy += abs2(p.field[j, i])
        end
    end
    return energy * kpix^2
end

# Check finite, nonzero samples without squaring values that can underflow or overflow.
function _pupil_field_scale(p::PupilFunction{T}) where {T}
    scale = zero(T)
    for i in axes(p.field, 2), j in axes(p.field, 1)
        _, _, _, inside = _pupil_sample(p, i, j)
        if inside
            value = p.field[j, i]
            scale = max(scale, abs(real(value)), abs(imag(value)))
        end
    end
    return scale
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

    for i in 1:sz, j in 1:sz
        _, _, kr2, inside = _pupil_sample(p, i, j)
        if inside
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
    all_float32 = T === Float32 && x isa Float32 && y isa Float32 && z isa Float32
    C = all_float32 ? Float32 : Float64
    R = promote_type(C, T, typeof(x), typeof(y), typeof(z))
    if isnan(x) || isnan(y) || isnan(z)
        return Complex{R}(NaN, NaN)
    end
    halfwidth = _alias_free_halfwidth(p)
    if _outside_field_window(x, halfwidth) || _outside_field_window(y, halfwidth)
        return zero(Complex{R})
    end
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    # Keep Float32 phase arithmetic, with enough range for the unweighted field sum.
    A = all_float32 ? Float64 : R
    result = zero(Complex{A})
    k₀² = k₀(p)^2

    @inbounds for i in 1:sz, j in 1:sz
        kx, ky, kr2, inside = _pupil_sample(p, i, j)
        if inside
            if all_float32
                kx, ky, kr2 = Float32(kx), Float32(ky), Float32(kr2)
            end
            delta = k₀² - kr2
            field = all_float32 ? ComplexF64(p.field[j, i]) : p.field[j, i]
            if all_float32 && delta >= 0
                kz = complex(sqrt(delta))
                phase = 2C(π) * (x * kx + y * ky + z * kz)
                result += field * exp(im * phase)
            else
                kz = sqrt(complex(delta))
                phase = 2C(π) * (x * kx + y * ky + z * kz)
                propagator = all_float32 ? exp(im * ComplexF64(phase)) : exp(im * phase)
                result += field * propagator
            end
        end
    end

    return all_float32 ? ComplexF32(result * Float64(kpix)^2) : result * kpix^2
end
