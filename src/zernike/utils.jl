# src/zernike/utils.jl

"""
Utility functions for working with Zernike polynomials and coefficients.
Includes common operations, aberration helpers, and manipulation functions.
"""

"""
    reset!(coeffs::ZernikeCoefficients) -> ZernikeCoefficients

Reset coefficients to default state (mag[1]=1, all others zero).
"""
function reset!(coeffs::ZernikeCoefficients)
    fill!(coeffs.mag, 0)
    fill!(coeffs.phase, 0)
    coeffs.mag[1] = 1
    return coeffs
end

"""
    add_aberration!(coeffs::ZernikeCoefficients, 
                   index::Integer, 
                   mag_value::Real=0, 
                   phase_value::Real=0;
                   indexing::ZernikeIndexing=OSA) -> ZernikeCoefficients

Add aberration terms to existing coefficients.

# Arguments
- `coeffs`: ZernikeCoefficients to modify
- `index`: Index of Zernike term (starting at 0)
- `mag_value`: Value to add to magnitude coefficient
- `phase_value`: Value to add to phase coefficient
- `indexing`: Indexing convention for the index

# Notes
- Index is 0-based to match standard Zernike notation
- Returns the modified coefficients for chaining
"""
function add_aberration!(coeffs::ZernikeCoefficients,
    index::Integer,
    mag_value::Real=0,
    phase_value::Real=0;
    indexing::ZernikeIndexing=OSA)
    index ≥ 0 || throw(ArgumentError("Index must be non-negative"))
    array_index = index + 1  # Convert to 1-based indexing for array access

    if array_index > length(coeffs.mag)
        throw(ArgumentError("Index $index exceeds coefficient array length"))
    end

    coeffs.mag[array_index] += mag_value
    coeffs.phase[array_index] += phase_value
    return coeffs
end

# Common aberration helpers
"""
    add_defocus!(coeffs::ZernikeCoefficients, amount::Real) -> ZernikeCoefficients

Add defocus aberration (n=2, l=0).
"""
function add_defocus!(coeffs::ZernikeCoefficients, amount::Real)
    return add_aberration!(coeffs, 4, 0, amount, indexing=OSA)
end

"""
    add_astigmatism!(coeffs::ZernikeCoefficients, 
                     amount::Real, 
                     angle::Real=0.0) -> ZernikeCoefficients

Add astigmatism with specified magnitude and angle.
"""
function add_astigmatism!(coeffs::ZernikeCoefficients, amount::Real, angle::Real=0.0)
    mag = amount * cos(2angle)
    add_aberration!(coeffs, 3, 0, mag, indexing=OSA)  # 0° astigmatism
    mag = amount * sin(2angle)
    add_aberration!(coeffs, 5, 0, mag, indexing=OSA)  # 45° astigmatism
    return coeffs
end

"""
    add_coma!(coeffs::ZernikeCoefficients, 
              amount::Real, 
              angle::Real=0.0) -> ZernikeCoefficients

Add coma with specified magnitude and angle.
"""
function add_coma!(coeffs::ZernikeCoefficients, amount::Real, angle::Real=0.0)
    mag = amount * cos(angle)
    add_aberration!(coeffs, 6, 0, mag, indexing=OSA)  # Vertical coma
    mag = amount * sin(angle)
    add_aberration!(coeffs, 7, 0, mag, indexing=OSA)  # Horizontal coma
    return coeffs
end

"""
    add_spherical!(coeffs::ZernikeCoefficients, amount::Real) -> ZernikeCoefficients

Add primary spherical aberration.
"""
function add_spherical!(coeffs::ZernikeCoefficients, amount::Real)
    return add_aberration!(coeffs, 11, 0, amount, indexing=OSA)
end

# Coefficient manipulation
"""
    scale!(coeffs::ZernikeCoefficients, factor::Real) -> ZernikeCoefficients

Scale all coefficients (except piston) by given factor.
"""
function scale!(coeffs::ZernikeCoefficients, factor::Real)
    for i in 2:length(coeffs.mag)
        coeffs.mag[i] *= factor
        coeffs.phase[i] *= factor
    end
    return coeffs
end

"""
    merge!(target::ZernikeCoefficients, 
           source::ZernikeCoefficients, 
           weight::Real=1.0) -> ZernikeCoefficients

Merge source coefficients into target with optional weighting.
"""
function merge!(target::ZernikeCoefficients,
    source::ZernikeCoefficients,
    weight::Real=1.0)
    n = min(length(target.mag), length(source.mag))
    for i in 1:n
        target.mag[i] += weight * source.mag[i]
        target.phase[i] += weight * source.phase[i]
    end
    return target
end

"""
    rms(coeffs::ZernikeCoefficients) -> Tuple{Float64,Float64}

Calculate RMS values for magnitude and phase coefficients (excluding piston).
"""
function rms(coeffs::ZernikeCoefficients)
    mag_rms = sqrt(sum(abs2, @view coeffs.mag[2:end]) / (length(coeffs.mag) - 1))
    phase_rms = sqrt(sum(abs2, @view coeffs.phase[2:end]) / (length(coeffs.phase) - 1))
    return (mag_rms, phase_rms)
end

"""
    trim!(coeffs::ZernikeCoefficients, threshold::Real) -> ZernikeCoefficients

Set coefficients below threshold (relative to RMS) to zero.
"""
function trim!(coeffs::ZernikeCoefficients, threshold::Real)
    mag_rms, phase_rms = rms(coeffs)

    for i in 2:length(coeffs.mag)
        if abs(coeffs.mag[i]) < threshold * mag_rms
            coeffs.mag[i] = 0
        end
        if abs(coeffs.phase[i]) < threshold * phase_rms
            coeffs.phase[i] = 0
        end
    end
    return coeffs
end

# Reporting
"""
    significant_terms(coeffs::ZernikeCoefficients, 
                     threshold::Real=0.01) -> Vector{Tuple{Int,Float64,Float64}}

Return list of significant terms: (index, magnitude, phase) above threshold.
"""
function significant_terms(coeffs::ZernikeCoefficients, threshold::Real=0.01)
    mag_rms, phase_rms = rms(coeffs)
    terms = Tuple{Int,Float64,Float64}[]

    for i in 1:length(coeffs.mag)
        if abs(coeffs.mag[i]) > threshold * mag_rms ||
           abs(coeffs.phase[i]) > threshold * phase_rms
            push!(terms, (i - 1, coeffs.mag[i], coeffs.phase[i]))
        end
    end

    return sort(terms, by=x -> max(abs(x[2]) / mag_rms, abs(x[3]) / phase_rms),
        rev=true)
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


