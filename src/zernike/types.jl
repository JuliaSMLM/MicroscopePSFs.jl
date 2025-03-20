# src/zernike/types.jl

"""
    ZernikeCoefficients{T<:Real}

Mutable structure to hold Zernike coefficients for both magnitude and phase of a pupil function.
Uses Noll indexing convention starting at index 1.

# Fields
- `mag::Vector{T}`: Coefficients for magnitude (1-indexed per Noll convention)
- `phase::Vector{T}`: Coefficients for phase (1-indexed per Noll convention)

# Notes
- First coefficient (index 1) typically represents piston
- Magnitude coefficients are typically normalized with mag[1] = 1
- Phase coefficients represent phase in radians
- RMS normalization (Noll convention) is used so coefficients directly represent RMS wavefront error
"""
mutable struct ZernikeCoefficients{T<:Real}
    mag::Vector{T}
    phase::Vector{T}
    
    function ZernikeCoefficients(mag::Vector{T}, phase::Vector{T}) where T<:Real
        length(mag) == length(phase) || 
            throw(ArgumentError("Magnitude and phase vectors must have same length"))
        new{T}(mag, phase)
    end
end

"""
    ZernikeCoefficients(n::Integer; T::Type{<:Real}=Float64)

Construct ZernikeCoefficients with `n` terms initialized to standard values.

# Arguments
- `n`: Number of Zernike terms
- `T`: Numeric type for coefficients (default: Float64)

# Returns
- `ZernikeCoefficients` with magnitude[1] = 1 and all other terms zero

# Examples
```julia
# Create coefficients for first 15 Zernike terms
zc = ZernikeCoefficients(15)

# Create coefficients with specific numeric type
zc = ZernikeCoefficients(10, T=Float32)
```
"""
function ZernikeCoefficients(n::Integer; T::Type{<:Real}=Float64)
    n > 0 || throw(ArgumentError("Number of coefficients must be positive"))
    mag = zeros(T, n)
    mag[1] = one(T)
    return ZernikeCoefficients(mag, zeros(T, n))
end

# Convenience methods for type inspection
Base.eltype(::ZernikeCoefficients{T}) where T = T
Base.length(zc::ZernikeCoefficients) = length(zc.mag)

# Indexing methods
Base.getindex(zc::ZernikeCoefficients, i::Integer) = (mag=zc.mag[i], phase=zc.phase[i])
Base.setindex!(zc::ZernikeCoefficients, v::Real, i::Integer, field::Symbol) = setfield!(zc, field)[i] = v

# Show methods for nice printing
function Base.show(io::IO, zc::ZernikeCoefficients{T}) where T
    print(io, "ZernikeCoefficients{$T}($(length(zc)) terms)")
end

function Base.show(io::IO, ::MIME"text/plain", zc::ZernikeCoefficients)
    println(io, "ZernikeCoefficients with $(length(zc)) terms (Noll indexed):")
    println(io, "  Magnitude coefficients:")
    for i in eachindex(zc.mag)
        abs(zc.mag[i]) > 1e-10 && println(io, "    [$i] = $(zc.mag[i])")
    end
    println(io, "  Phase coefficients:")
    for i in eachindex(zc.phase)
        abs(zc.phase[i]) > 1e-10 && println(io, "    [$i] = $(zc.phase[i])")
    end
end

# Iteration interface
Base.iterate(zc::ZernikeCoefficients) = iterate((zc.mag, zc.phase))
Base.iterate(zc::ZernikeCoefficients, state) = iterate((zc.mag, zc.phase), state)

# Equality comparison
function Base.:(==)(a::ZernikeCoefficients, b::ZernikeCoefficients)
    return a.mag == b.mag && a.phase == b.phase
end

# Copying
Base.copy(zc::ZernikeCoefficients) = ZernikeCoefficients(copy(zc.mag), copy(zc.phase))