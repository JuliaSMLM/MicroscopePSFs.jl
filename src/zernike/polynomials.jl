# src/zernike/polynomials.jl

"""
    radialpolynomial(n::Integer, m::Integer, ρ::Real) -> Real

Compute the unnormalized radial component R_n^m(ρ) of the Zernike polynomial.

# Arguments
- `n`: Radial order
- `m`: Azimuthal order (absolute value of azimuthal frequency)
- `ρ`: Radial coordinate (0 ≤ ρ ≤ 1)

# Notes
- Returns 0 for ρ > 1
- No normalization applied - returns the standard mathematical form
"""
function radialpolynomial(n::Integer, m::Integer, ρ::Real)
    # Input validation
    n ≥ 0 || throw(ArgumentError("Radial order n must be non-negative"))
    m ≥ 0 || throw(ArgumentError("Azimuthal order m must be non-negative"))
    m ≤ n || throw(ArgumentError("Azimuthal order m must be ≤ n"))
    mod(n - m, 2) == 0 || throw(ArgumentError("n - m must be even"))
    
    # Early return for points outside unit circle
    ρ > 1 && return zero(ρ)
    
    result = zero(typeof(ρ))
    
    # Sum over k
    for k in 0:((n-m)÷2)
        # Correctly calculate the factorial (n-k)!
        num = (-1)^k
        for i in 1:(n-k)
            num *= i
        end
        
        den = factorial(k) * factorial((n+m)÷2 - k) * factorial((n-m)÷2 - k)
        result += (num/den) * ρ^(n-2k)
    end
    
    return result
end

"""
    zernikepolynomial(n::Integer, l::Integer, ρ::Real, ϕ::Real) -> Real

Compute the complete Zernike polynomial Z_n^l(ρ,ϕ) with Noll normalization (RMS=1).

# Arguments
- `n`: Radial order
- `l`: Azimuthal frequency (signed)
- `ρ`: Radial coordinate (0 ≤ ρ ≤ 1)
- `ϕ`: Azimuthal angle in radians

# Notes
- Uses Noll normalization where RMS=1 over unit circle
- Combines radial polynomial with appropriate trigonometric function
"""
function zernikepolynomial(n::Integer, l::Integer, ρ::Real, ϕ::Real)
    m = abs(l)
    
    # Normalization factor for Noll (RMS = 1)
    norm_factor = sqrt(2 * (n + 1))
    if m == 0
        norm_factor = sqrt(n + 1)
    end
    
    # Get standard radial polynomial
    r = radialpolynomial(n, m, ρ)
    
    # Apply normalization and angular component
    if l < 0
        return norm_factor * r * sin(m*ϕ)
    else
        return norm_factor * r * cos(m*ϕ)
    end
end

"""
    zernikepolynomial(j::Integer, ρ::Real, ϕ::Real) -> Real

Compute Zernike polynomial using Noll index.

# Arguments
- `j`: Polynomial index using Noll convention
- `ρ`: Radial coordinate
- `ϕ`: Azimuthal angle
"""
function zernikepolynomial(j::Integer, ρ::Real, ϕ::Real)
    n, l = noll2nl(j)
    return zernikepolynomial(n, l, ρ, ϕ)
end

"""
    evaluate_pupil(coeffs::ZernikeCoefficients, grid_size::Integer) -> Matrix{Complex{Float64}}

Generate complex pupil function from Zernike coefficients.

# Arguments
- `coeffs`: ZernikeCoefficients containing magnitude and phase coefficients
- `grid_size`: Size of the output grid (grid_size × grid_size)

# Returns
- Complex-valued matrix representing the pupil function

# Notes
- Output grid is normalized to unit circle
- Points outside unit circle are set to zero
- Phase is applied as exp(iϕ)
- Uses Noll indexing throughout
"""
function evaluate_pupil(coeffs::ZernikeCoefficients, grid_size::Integer)
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
        for k in eachindex(coeffs.mag)
            basis = zernikepolynomial(k, ρ, ϕ)
            mag += coeffs.mag[k] * basis
            phase += coeffs.phase[k] * basis
        end
        
        pupil[j,i] = mag * exp(im * phase)
    end
    
    return pupil
end

"""
    max_radial_order(num_coeffs::Integer)

Calculate maximum radial order N given number of coefficients L.
Solves (N+1)(N+2)/2 = L
"""
function max_radial_order(num_coeffs::Integer)
    # Solve quadratic equation (N+1)(N+2)/2 = L
    # N² + 3N + (2-2L) = 0
    a = 1
    b = 3
    c = 2 - 2 * num_coeffs
    N = (-b + sqrt(b^2 - 4a * c)) / (2a)
    return floor(Int, N)
end