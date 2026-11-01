# src/emitters.jl

import SMLMData: AbstractEmitter 

"""
    DipoleVector{T<:Real}

A normalized 3D dipole orientation. Every constructor normalizes nonzero inputs.
All input components must be finite. Normalization uses a working type promoted
from the floating-point types of every input and the requested component type.
Conversion to the requested type must preserve a finite, nonzero direction.
The zero vector is allowed for rotating-PSF metadata; fixed dipoles must be nonzero.
The untyped factory promotes integers and rationals to floating point.
`DipoleVector{T}(px, py, pz)` retains `T`. Integer and rational types require
exactly representable normalized components; otherwise use a floating-point type
or the untyped `DipoleVector(px, py, pz)`. Integer types allow unit axes and zero.

# Fields
- `px::T`: x component of the dipole vector
- `py::T`: y component of the dipole vector
- `pz::T`: z component of the dipole vector
"""
struct DipoleVector{T<:Real}
    px::T
    py::T
    pz::T

    function DipoleVector{T}(px::Real, py::Real, pz::Real) where {T<:Real}
        all(isfinite, (px, py, pz)) || throw(ArgumentError(
            "Dipole components must be finite; supply finite orientation components"))
        if T <: Union{Integer,Rational}
            x, y, z = _exact_dipole_components(T, px, py, pz)
            return new{T}(x, y, z)
        end
        W = promote_type(float(typeof(px)), float(typeof(py)), float(typeof(pz)), float(T))
        x, y, z = W(px), W(py), W(pz)
        scale = max(abs(x), abs(y), abs(z))
        iszero(scale) && return new{T}(T(x), T(y), T(z))
        x, y, z = x / scale, y / scale, z / scale
        magnitude = sqrt(x^2 + y^2 + z^2)
        normalized = T.((x / magnitude, y / magnitude, z / magnitude))
        all(isfinite, normalized) && !all(iszero, normalized) || throw(ArgumentError(
            "Normalized dipole conversion must preserve a finite nonzero direction; " *
            "use a wider floating-point type or the untyped DipoleVector(px, py, pz)"))
        return new{T}(normalized...)
    end
end

function DipoleVector(px::Real, py::Real, pz::Real)
    T = float(promote_type(typeof(px), typeof(py), typeof(pz)))
    return DipoleVector{T}(px, py, pz)
end

function _dipole_type_error(T)
    return ArgumentError(
        "Normalized dipole components are not exactly representable in $T; " *
        "use a floating-point type or the untyped DipoleVector(px, py, pz)")
end

# Exact arithmetic distinguishes rational unit vectors from irrational directions.
function _exact_dipole_components(::Type{T}, px, py, pz) where {T}
    q = Rational{BigInt}.((px, py, pz))
    norm² = sum(abs2, q)
    iszero(norm²) && return (zero(T), zero(T), zero(T))
    a, b = numerator(norm²), denominator(norm²)
    na, nb = isqrt(a), isqrt(b)
    na^2 == a && nb^2 == b || throw(_dipole_type_error(T))
    normalized = q ./ (na // nb)
    try
        result = T.(normalized)
        result == normalized || throw(_dipole_type_error(T))
        return result
    catch e
        e isa InexactError || e isa OverflowError || rethrow()
        throw(_dipole_type_error(T))
    end
end

function _require_fixed_dipole(dipole::DipoleVector)
    all(iszero, (dipole.px, dipole.py, dipole.pz)) && throw(ArgumentError(
        "A fixed dipole must be nonzero; supply a nonzero dipole orientation"))
    return dipole
end

"""
    DipoleEmitter3D{T} <: AbstractEmitter

3D dipole emitter with position, orientation and optical properties.

`DipoleEmitter3D{T}(x, y, z, photons, dipole)` retains `T` when its normalized
orientation and scalar arguments are representable in `T`. Otherwise use a
floating-point type or the untyped `DipoleEmitter3D(x, y, z, photons, dipole)`,
which promotes integer and rational inputs to floating point.

# Fields
- `x::T`: x-coordinate in microns
- `y::T`: y-coordinate in microns 
- `z::T`: z-coordinate in microns
- `photons::T`: number of photons
- `dipole::DipoleVector{T}`: dipole orientation vector
"""
mutable struct DipoleEmitter3D{T} <: AbstractEmitter
    x::T
    y::T
    z::T 
    photons::T
    dipole::DipoleVector{T}

    function DipoleEmitter3D{T}(x::Real, y::Real, z::Real, photons::Real,
        dipole::DipoleVector) where {T<:Real}
        _require_fixed_dipole(dipole)
        try
            orientation = DipoleVector{T}(dipole.px, dipole.py, dipole.pz)
            return new{T}(T(x), T(y), T(z), T(photons), orientation)
        catch e
            if e isa InexactError || e isa OverflowError || e isa ArgumentError
                throw(ArgumentError(
                    "Emitter values or normalized dipole cannot be represented in $T; " *
                    "use a floating-point type or the untyped " *
                    "DipoleEmitter3D(x, y, z, photons, dipole)"))
            end
            rethrow()
        end
    end
end

"""
    DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real, dipole::DipoleVector)

Create an emitter at `(x, y, z)` in microns with the given photon count and a
normalized, nonzero dipole orientation. Numeric arguments promote together;
integer and rational inputs promote to floating point.
Returns a `DipoleEmitter3D`; use a nonzero `dipole` for a fixed emitter.

# Examples
```julia
emitter = DipoleEmitter3D(1.0, 1.0, 0.0, 1000.0, DipoleVector(2.0, 0.0, 0.0))
```
"""
function DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real,
    dipole::DipoleVector)
    T = float(promote_type(typeof(x), typeof(y), typeof(z), typeof(photons),
        typeof(dipole.px)))
    return DipoleEmitter3D{T}(x, y, z, photons, dipole)
end

# Constructor with vector normalization
"""
    DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real,
                    dx::Real, dy::Real, dz::Real)

Create a `DipoleEmitter3D` with specified position, number of photons, and dipole orientation.
The dipole orientation vector will be normalized and must be nonzero.

# Arguments
- `x::Real`: x-coordinate in microns
- `y::Real`: y-coordinate in microns
- `z::Real`: z-coordinate in microns
- `photons::Real`: number of photons
- `dx::Real`: x component of dipole orientation (normalized)
- `dy::Real`: y component of dipole orientation (normalized)
- `dz::Real`: z component of dipole orientation (normalized)
"""
function DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real,
                        dx::Real, dy::Real, dz::Real)
    return DipoleEmitter3D(x, y, z, photons, DipoleVector(dx, dy, dz))
end



