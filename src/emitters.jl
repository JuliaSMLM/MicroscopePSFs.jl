import SMLMData: AbstractEmitter 

"""
    DipoleVector{T} <: Real

A 3D vector representing the dipole orientation.

# Fields
- `px::T`: x component of the dipole vector
- `py::T`: y component of the dipole vector
- `pz::T`: z component of the dipole vector
"""
struct DipoleVector{T<:Real}
    px::T
    py::T
    pz::T
end

function DipoleVector(px::Real, py::Real, pz::Real)
    T = promote_type(typeof.((px, py, pz))...)
    
    # Normalize dipole vector
    d = [px, py, pz]
    d_norm = norm(d)
    d_norm > 0 || throw(ArgumentError("dipole vector cannot be zero"))
    d = d ./ d_norm
    
    DipoleVector{T}(T(d[1]), T(d[2]), T(d[3]))
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
end

# Constructor with vector normalization
"""
    DipoleEmitter3D(x::Real, y::Real, z::Real, photons::Real,
                    dx::Real, dy::Real, dz::Real)

Create a `DipoleEmitter3D` with specified position, number of photons, and dipole orientation.
The dipole orientation vector will be normalized.

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
    T = promote_type(typeof.((x,y,z,photons,dx,dy,dz))...)
    
    dipole = DipoleVector(dx, dy, dz)
    
    DipoleEmitter3D{T}(T(x), T(y), T(z), T(photons), dipole)
end



