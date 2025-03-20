# src/zernike/utils.jl

"""
Utility functions for working with Zernike polynomials and coefficients.
Provides core analysis functions.
Uses Noll indexing throughout.
"""

"""
    rms(coeffs::ZernikeCoefficients) -> Tuple{Float64,Float64}

Calculate RMS values for magnitude and phase coefficients (excluding piston).
"""
function rms(coeffs::ZernikeCoefficients)
    # Get indices excluding piston (index 1)
    indices = [i for i in eachindex(coeffs.mag) if i > 1]
    
    # Calculate RMS excluding piston
    mag_rms = sqrt(sum(abs2, coeffs.mag[indices]) / length(indices))
    phase_rms = sqrt(sum(abs2, coeffs.phase[indices]) / length(indices))
    
    return (mag_rms, phase_rms)
end

"""
    significant_terms(coeffs::ZernikeCoefficients, 
                     threshold::Real=0.01) -> Vector{Tuple{Int,Float64,Float64}}

Return list of significant terms: (index, magnitude, phase) above threshold.
"""
function significant_terms(coeffs::ZernikeCoefficients, threshold::Real=0.01)
    mag_rms, phase_rms = rms(coeffs)
    terms = Tuple{Int,Float64,Float64}[]

    for i in eachindex(coeffs.mag)
        if abs(coeffs.mag[i]) > threshold * mag_rms ||
           abs(coeffs.phase[i]) > threshold * phase_rms
            push!(terms, (i, coeffs.mag[i], coeffs.phase[i]))
        end
    end

    return sort(terms, by=x -> max(abs(x[2]) / mag_rms, abs(x[3]) / phase_rms),
        rev=true)
end