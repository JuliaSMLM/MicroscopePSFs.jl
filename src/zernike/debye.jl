# src/zernike/debye.jl

"""
    propagate_field(coeffs::ZernikeCoefficients{T},
                   λ::T, n::T, nₐ::T,
                   r::T, φ::T, z::T;
                   nk::Integer=100) where {T<:Real}

Propagate field from pupil described by Zernike coefficients to point (r,φ,z).

U(r,φ,z) = 2π ∫₀^k_max dk_r k_r e^{2πiz√(k²-k_r²)} ∑_{n,m} R_n^m J_m(2πk_r r)[C_a cos(mφ) + C_b sin(mφ)]

C_a is the complex coefficient for the Zernike mode (n,m) and C_b is the complex coefficient for (n,-m).
"""
function propagate_field(coeffs::ZernikeCoefficients{T},
    λ::T, n::T, nₐ::T,
    r::T, φ::T, z::T;
    nk::Integer=100) where {T<:Real}

    # Wave parameters
    k₀ = 2π * n / λ
    k_max = 2π * nₐ / λ

    # Integration grid for k_r
    k_r = range(zero(T), stop=k_max, length=nk)
    dk = k_max / (nk - 1)

    # Initialize field
    field = zero(Complex{T})

    # Maximum order based on number of coefficients
    max_n = max_radial_order(length(coeffs.mag))

    # Loop over k_r integration points
    for (i, kr) in enumerate(k_r)
        kr > k_max && continue

        # Calculate kz component
        kz = sqrt(complex(k₀^2 - kr^2))

        # Get radial coordinate in pupil (normalized)
        ρ = kr / k_max

        # Initialize pupil function at this kr
        pupil = zero(Complex{T})

        # Sum over Zernike orders
        for n in 0:max_n
            for m in 0:n
                # Skip if n-m is odd
                (n - m) % 2 == 1 && continue

                # Get radial polynomial
                R = radialpolynomial(n, m, ρ)

                # Calculate coefficient indices in Noll ordering
                a = nl2osa(n, m) + 1 # +1 for 1-based indexing
                b = nl2osa(n, -m) + 1

                # Skip if beyond coefficient array
                (b ≥ length(coeffs)) && continue

                # print if coefficients are non zero
                if coeffs.mag[a] != 0 || coeffs.mag[b] != 0
                    println("n: $n, m: $m, a: $a, b: $b, mag: $(coeffs.mag[a]), $(coeffs.mag[b])")
                end

                if coeffs.phase[a] != 0 || coeffs.phase[b] != 0
                    println("n: $n, m: $m, a: $a, b: $b, phase: $(coeffs.phase[a]), $(coeffs.phase[b])")
                end


                # Get coefficients in U + iV form
                Ua = coeffs.mag[a] * cos(coeffs.phase[a])
                Va = coeffs.mag[a] * sin(coeffs.phase[a])
                Ca = Complex(Ua, Va)

                Ub = coeffs.mag[b] * cos(coeffs.phase[b])
                Vb = coeffs.mag[b] * sin(coeffs.phase[b])
                Cb = Complex(Ub, Vb)

                # Get Bessel function of order m
                bessel_term = besselj(m, 2π * kr * r)

                if m == 0
                    pupil += R * bessel_term * Ca
                else
                    pupil += R * bessel_term * (Ca * cos(m * φ) + Cb * sin(m * φ))
                end
            end
        end

        # Add propagation phase
        pupil *= exp(2π * im * z * kz)

        # Add to integral with kr factor and weights
        weight = (i == 1 || i == nk) ? 0.5 : 1.0
        field += weight * kr * pupil * dk
    end

    # Apply prefactor
    field *= 2π

    return field
end