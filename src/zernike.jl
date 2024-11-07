# Define indexing convention as enum instead of strings
@enum ZernikeIndexing OSA Noll

"""
    ZernikeCoefficients{T<:Real}

Mutable structure to hold Zernike coefficients for magnitude and phase.

# Fields
- `mag`: Zernike coefficients of pupil magnitude
- `phase`: Zernike coefficients of pupil phase
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

# Convenience constructor
function ZernikeCoefficients(n::Integer; T::Type{<:Real}=Float64)
    mag = zeros(T, n)
    mag[1] = one(T)
    return ZernikeCoefficients(mag, zeros(T, n))
end

# Helper methods
function reset!(coeffs::ZernikeCoefficients)
    fill!(coeffs.mag, 0)
    fill!(coeffs.phase, 0)
    coeffs.mag[1] = 1
    return coeffs
end

function add_aberration!(coeffs::ZernikeCoefficients, index::Integer, 
                        mag_value::Real=0, phase_value::Real=0)
    coeffs.mag[index] += mag_value
    coeffs.phase[index] += phase_value
    return coeffs
end

# Index conversion functions with enum
"""
    zernikepolynomial(j::Integer, ρ::Real, ϕ::Real, indexing::ZernikeIndexing=OSA)

Compute jth Zernike polynomial using specified indexing convention.
"""
function zernikepolynomial(j::Integer, ρ::Real, ϕ::Real, indexing::ZernikeIndexing=OSA)
    n, l = if indexing == OSA
        osa2nl(j)
    else
        noll2nl(j)
    end
    return zernikepolynomial(n, l, ρ, ϕ)
end

# Main computational functions without caching
function radialpolynomial(n::Integer, m::Integer, ρ::Real)
    ρ > 1 && return zero(ρ)
    
    g = m == 0 ? sqrt(n + 1) : sqrt(2n + 2)
    
    result = zero(promote_type(typeof(g), typeof(ρ)))
    for k in 0:div(n-m, 2)
        coef = g * (-1)^k * 
               prod((n - m)÷2 - k + 1 : n - k) / 
               (factorial(k) * factorial((n + m)÷2 - k))
        result += coef * ρ^(n-2k)
    end
    
    return result
end

function zernikepolynomial(n::Integer, l::Integer, ρ::Real, ϕ::Real)
    m = abs(l)
    r = radialpolynomial(n, m, ρ)
    return l < 0 ? r * sin(m*ϕ) : r * cos(m*ϕ)
end

# Pupil evaluation function
"""
    evaluate_pupil(coeffs::ZernikeCoefficients, grid_size::Integer; 
                  indexing::ZernikeIndexing=OSA)

Generate complex pupil function from Zernike coefficients.
"""
function evaluate_pupil(coeffs::ZernikeCoefficients, grid_size::Integer;
                       indexing::ZernikeIndexing=OSA)
    T = eltype(coeffs.mag)
    pupil = zeros(Complex{T}, grid_size, grid_size)
    
    xs = ys = range(-1, 1, length=grid_size)
    
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        ϕ = atan(y, x)
        
        mag = phase = zero(T)
        for k in 1:length(coeffs.mag)
            basis = zernikepolynomial(k-1, ρ, ϕ, indexing)
            mag += coeffs.mag[k] * basis
            phase += coeffs.phase[k] * basis
        end
        
        pupil[j,i] = mag * exp(im * phase)
    end
    
    return pupil
end