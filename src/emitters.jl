# src/emitters.jl

import SMLMData: AbstractEmitter 

"""
    DipoleVector{T<:Real}

A normalized 3D dipole orientation. Every constructor normalizes nonzero inputs.
All components must be finite and remain finite after conversion to the requested type.
The zero vector is allowed for rotating-PSF metadata; fixed dipoles must be nonzero.
Integer and rational component types promote to floating point for normalization.

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
        x, y, z = float(convert(T, px)), float(convert(T, py)), float(convert(T, pz))
        all(isfinite, (x, y, z)) || throw(ArgumentError(
            "Dipole components must remain finite after conversion; " *
            "use representable finite components"))
        scale = max(abs(x), abs(y), abs(z))
        iszero(scale) && return new{typeof(x)}(x, y, z)
        x, y, z = x / scale, y / scale, z / scale
        magnitude = sqrt(x^2 + y^2 + z^2)
        return new{typeof(x / magnitude)}(x / magnitude, y / magnitude, z / magnitude)
    end
end

function DipoleVector(px::Real, py::Real, pz::Real)
    T = promote_type(typeof(px), typeof(py), typeof(pz))
    return DipoleVector{T}(px, py, pz)
end

function _require_fixed_dipole(dipole::DipoleVector)
    all(iszero, (dipole.px, dipole.py, dipole.pz)) && throw(ArgumentError(
        "A fixed dipole must be nonzero; supply a nonzero dipole orientation"))
    return dipole
end

"""
    DipoleEmitter3D{T} <: AbstractEmitter

3D dipole emitter with position, orientation and optical properties.

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
        orientation = _require_fixed_dipole(
            DipoleVector{T}(dipole.px, dipole.py, dipole.pz))
        U = typeof(orientation.px)
        return new{U}(U(x), U(y), U(z), U(photons), orientation)
    end
end

"""
    DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real, dipole::DipoleVector)

Create an emitter at `(x, y, z)` in microns with the given photon count and a
normalized, nonzero dipole orientation. Numeric arguments promote together.
Returns a `DipoleEmitter3D`; use a nonzero `dipole` for a fixed emitter.

# Examples
```julia
emitter = DipoleEmitter3D(1.0, 1.0, 0.0, 1000.0, DipoleVector(2.0, 0.0, 0.0))
```
"""
function DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real,
    dipole::DipoleVector)
    T = promote_type(typeof(x), typeof(y), typeof(z), typeof(photons), typeof(dipole.px))
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



