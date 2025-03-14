# src/zernike/polynomials.jl

"""
Implementations of Zernike polynomial calculations and pupil function generation.
Includes radial polynomials, full Zernike polynomials, and pupil function evaluation.
"""

"""
    radialpolynomial(n::Integer, m::Integer, ρ::Real) -> Real

Compute the radial component R_n^m(ρ) of the Zernike polynomial.

# Arguments
- `n`: Radial order
- `m`: Azimuthal order (absolute value of azimuthal frequency)
- `ρ`: Radial coordinate (0 ≤ ρ ≤ 1)

# Notes
- Returns 0 for ρ > 1
- Normalizes according to born-wolf convention
"""
function radialpolynomial(n::Integer, m::Integer, ρ::Real)
    # Input validation
    n ≥ 0 || throw(ArgumentError("Radial order n must be non-negative"))
    m ≥ 0 || throw(ArgumentError("Azimuthal order m must be non-negative"))
    m ≤ n || throw(ArgumentError("Azimuthal order m must be ≤ n"))
    mod(n - m, 2) == 0 || throw(ArgumentError("n - m must be even"))
    
    # Early return for points outside unit circle
    ρ > 1 && return zero(ρ)
    
    # Normalization factor
    g = m == 0 ? sqrt(n + 1) : sqrt(2n + 2)
    
    result = zero(promote_type(typeof(g), typeof(ρ)))
    
    # Sum over k
    for k in 0:((n-m)÷2)
        num = (-1)^k * prod((n-k):-1:(n-k-(k-1)))
        den = factorial(k) * factorial((n+m)÷2 - k) * factorial((n-m)÷2 - k)
        result += (num/den) * ρ^(n-2k)
    end
    
    return g * result
end

"""
    zernikepolynomial(n::Integer, l::Integer, ρ::Real, ϕ::Real) -> Real

Compute the complete Zernike polynomial Z_n^l(ρ,ϕ).

# Arguments
- `n`: Radial order
- `l`: Azimuthal frequency (signed)
- `ρ`: Radial coordinate (0 ≤ ρ ≤ 1)
- `ϕ`: Azimuthal angle in radians

# Notes
- Uses born-wolf normalization
- Combines radial polynomial with appropriate trigonometric function
"""
function zernikepolynomial(n::Integer, l::Integer, ρ::Real, ϕ::Real)
    m = abs(l)
    r = radialpolynomial(n, m, ρ)
    if l < 0
        return r * sin(m*ϕ)
    else
        return r * cos(m*ϕ)
    end
end

"""
    zernikepolynomial(j::Integer, ρ::Real, ϕ::Real, indexing::ZernikeIndexing=OSA) -> Real

Compute Zernike polynomial using single-index notation.

# Arguments
- `j`: Polynomial index (OSA or Noll)
- `ρ`: Radial coordinate
- `ϕ`: Azimuthal angle
- `indexing`: Index convention to use (OSA or Noll)
"""
function zernikepolynomial(j::Integer, ρ::Real, ϕ::Real, indexing::ZernikeIndexing=OSA)
    n, l = get_nl(j, indexing)
    return zernikepolynomial(n, l, ρ, ϕ)
end

"""
    evaluate_pupil(coeffs::ZernikeCoefficients, grid_size::Integer;
                  indexing::ZernikeIndexing=OSA) -> Matrix{Complex{Float64}}

Generate complex pupil function from Zernike coefficients.

# Arguments
- `coeffs`: ZernikeCoefficients containing magnitude and phase coefficients
- `grid_size`: Size of the output grid (grid_size × grid_size)
- `indexing`: Indexing convention for coefficients (OSA or Noll)

# Returns
- Complex-valued matrix representing the pupil function

# Notes
- Output grid is normalized to unit circle
- Points outside unit circle are set to zero
- Phase is applied as exp(iϕ)
"""
function evaluate_pupil(coeffs::ZernikeCoefficients, grid_size::Integer;
                       indexing::ZernikeIndexing=OSA)
    T = eltype(coeffs.mag)
    pupil = zeros(Complex{T}, grid_size, grid_size)
    
    # Create normalized coordinate grid
    xs = ys = range(-1, 1, length=grid_size)
    
    # Evaluate polynomials at each point
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        
        ϕ = atan(y, x)
        
        # Sum magnitude and phase contributions
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

"""
    evaluate_polynomial_grid(n::Integer, l::Integer, grid_size::Integer) -> Matrix{Float64}

Generate a grid of Zernike polynomial values for visualization.

# Arguments
- `n`: Radial order
- `l`: Azimuthal frequency
- `grid_size`: Size of output grid

# Returns
- Matrix of polynomial values over [-1,1] × [-1,1]
"""
function evaluate_polynomial_grid(n::Integer, l::Integer, grid_size::Integer)
    result = zeros(grid_size, grid_size)
    xs = ys = range(-1, 1, length=grid_size)
    
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ϕ = atan(y, x)
        result[j,i] = zernikepolynomial(n, l, ρ, ϕ)
    end
    
    return result
end
