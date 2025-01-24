# src/zernike/fitting.jl

"""
Functions for fitting complex pupil functions to Zernike polynomial bases.
Includes both magnitude and phase fitting capabilities.
"""

"""
    build_zernike_basis(grid_size::Integer, num_terms::Integer;
                       indexing::ZernikeIndexing=OSA) -> Matrix{Float64}

Build a matrix of Zernike polynomials for least-squares fitting.

# Arguments
- `grid_size`: Size of the sampling grid
- `num_terms`: Number of Zernike terms to include
- `indexing`: Indexing convention for polynomials

# Returns
- Matrix of size (valid_points, num_terms) containing Zernike values
"""
function build_zernike_basis(grid_size::Integer, num_terms::Integer;
                           indexing::ZernikeIndexing=OSA)
    # Create coordinate grid
    xs = ys = range(-1, 1, length=grid_size)
    
    # Count valid points (inside unit circle)
    valid_points = count(x^2 + y^2 ≤ 1 for x in xs, y in ys)
    
    # Allocate basis matrix
    basis = zeros(valid_points, num_terms)
    
    # Fill basis matrix
    idx = 1
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        
        ϕ = atan(y, x)
        
        for k in 0:(num_terms-1)
            basis[idx, k+1] = zernikepolynomial(k, ρ, ϕ, indexing)
        end
        idx += 1
    end
    
    return basis
end

"""
    fit_pupil_to_zernike(pupil::AbstractMatrix{<:Complex}, num_terms::Integer;
                        indexing::ZernikeIndexing=OSA,
                        mask::Union{Nothing,AbstractMatrix{Bool}}=nothing) -> ZernikeCoefficients

Fit a complex pupil function to Zernike polynomials.

# Arguments
- `pupil`: Complex-valued pupil function matrix
- `num_terms`: Number of Zernike terms to fit
- `indexing`: Indexing convention to use
- `mask`: Optional boolean mask for valid pupil points

# Returns
- ZernikeCoefficients containing fitted magnitude and phase coefficients

# Notes
- Fits magnitude and phase separately using least squares
- Handles phase wrapping in the fit
- Optional mask allows fitting specific regions
"""
function fit_pupil_to_zernike(pupil::AbstractMatrix{<:Complex}, num_terms::Integer;
                            indexing::ZernikeIndexing=OSA,
                            mask::Union{Nothing,AbstractMatrix{Bool}}=nothing)
    grid_size = size(pupil, 1)
    size(pupil, 1) == size(pupil, 2) || 
        throw(ArgumentError("Pupil must be square"))
    
    if !isnothing(mask)
        size(mask) == size(pupil) ||
            throw(ArgumentError("Mask size must match pupil size"))
    end
    
    # Build Zernike basis
    basis = build_zernike_basis(grid_size, num_terms; indexing=indexing)
    
    # Extract valid pupil points
    xs = ys = range(-1, 1, length=grid_size)
    valid_points = Int[]
    pupil_mag = Float64[]
    pupil_phase = Float64[]
    
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        
        if !isnothing(mask) && !mask[j,i]
            continue
        end
        
        push!(valid_points, length(pupil_mag) + 1)
        push!(pupil_mag, abs(pupil[j,i]))
        push!(pupil_phase, angle(pupil[j,i]))
    end
    
    # Fit magnitude
    mag_coeffs = basis[valid_points,:] \ pupil_mag
    
    # Unwrap phase and fit
    phase_unwrapped = unwrap_phase(pupil_phase)
    phase_coeffs = basis[valid_points,:] \ phase_unwrapped
    
    return ZernikeCoefficients(mag_coeffs, phase_coeffs)
end

"""
    unwrap_phase(phase::AbstractVector{<:Real}) -> Vector{Float64}

Unwrap phase values to handle 2π discontinuities.
"""
function unwrap_phase(phase::AbstractVector{<:Real})
    result = copy(phase)
    for i in 2:length(result)
        while result[i] - result[i-1] > π
            result[i] -= 2π
        end
        while result[i] - result[i-1] < -π
            result[i] += 2π
        end
    end
    return result
end

"""
    fit_error(pupil::AbstractMatrix{<:Complex}, coeffs::ZernikeCoefficients;
             indexing::ZernikeIndexing=OSA) -> Float64

Calculate RMS error between pupil and its Zernike fit.

# Returns
- RMS error normalized to mean pupil magnitude
"""
function fit_error(pupil::AbstractMatrix{<:Complex}, coeffs::ZernikeCoefficients;
                  indexing::ZernikeIndexing=OSA)
    # Regenerate fitted pupil
    fitted = evaluate_pupil(coeffs, size(pupil, 1); indexing=indexing)
    
    # Calculate error only within unit circle
    grid_size = size(pupil, 1)
    xs = ys = range(-1, 1, length=grid_size)
    
    error_sum = 0.0
    count = 0
    mean_mag = 0.0
    
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        
        error_sum += abs2(pupil[j,i] - fitted[j,i])
        mean_mag += abs(pupil[j,i])
        count += 1
    end
    
    return sqrt(error_sum / count) / (mean_mag / count)
end
