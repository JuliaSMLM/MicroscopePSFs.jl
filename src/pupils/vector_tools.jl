# src/vector_tools.jl

"""
Helper functions for vectorial PSF calculations including multi-layer interfaces,
Fresnel coefficients, and proper dipole field computations.
"""

"""
    calculate_wave_vectors(kr2::Real, λ::Real, 
                         n_medium::Real, n_coverslip::Real, n_immersion::Real)

Calculate z-components of wave vectors in all three media.

# Arguments
- `kr2`: Squared lateral spatial frequency
- `λ`: Wavelength in microns
- `n_medium`: Refractive index of the sample medium
- `n_coverslip`: Refractive index of the coverslip
- `n_immersion`: Refractive index of the immersion medium

# Returns
- Tuple of (kz_medium, kz_coverslip, kz_immersion): z-components of wave vectors in each medium
"""
function calculate_wave_vectors(kr2::Real, λ::Real,
    n_medium::Real, n_coverslip::Real, n_immersion::Real)
    
    k₀ = 1 / λ  

    # Adjust formulas to maintain the same physical values
    kz_medium = k₀ * sqrt(complex(n_medium^2 - kr2 * λ^2))
    kz_coverslip = k₀ * sqrt(complex(n_coverslip^2 - kr2 * λ^2))
    kz_immersion = k₀ * sqrt(complex(n_immersion^2 - kr2 * λ^2))

    return kz_medium, kz_coverslip, kz_immersion
end

"""
    calculate_interface_fresnel(kr2::Real, λ::Real, 
                              n1::Real, n2::Real)

Calculate Fresnel transmission coefficients for a single interface.

# Arguments
- `kr2`: Squared lateral spatial frequency
- `λ`: Wavelength in microns
- `n1`: Refractive index of first medium
- `n2`: Refractive index of second medium

# Returns
- Tuple of (Tp, Ts): p and s polarization transmission coefficients
"""
function calculate_interface_fresnel(kr2::Real, λ::Real, n1::Real, n2::Real)
    # Calculate cosines of angles in both media
    cosθ1 = sqrt(complex(1 - kr2 * λ^2 / n1^2))
    cosθ2 = sqrt(complex(1 - kr2 * λ^2 / n2^2))

    # Calculate transmission coefficients
    Tp = 2 * n1 * cosθ1 / (n1 * cosθ2 + n2 * cosθ1)
    Ts = 2 * n1 * cosθ1 / (n1 * cosθ1 + n2 * cosθ2)

    return Tp, Ts
end

"""
    calculate_fresnel_coefficients(kr2::Real, λ::Real, 
                                 n_medium::Real, n_coverslip::Real, n_immersion::Real)

Calculate combined Fresnel transmission coefficients across all interfaces.

# Arguments
- `kr2`: Squared lateral spatial frequency
- `λ`: Wavelength in microns
- `n_medium`: Refractive index of the sample medium
- `n_coverslip`: Refractive index of the coverslip
- `n_immersion`: Refractive index of the immersion medium

# Returns
- Tuple of (Tp, Ts): Combined p and s polarization transmission coefficients
"""
function calculate_fresnel_coefficients(kr2::Real, λ::Real,
    n_medium::Real, n_coverslip::Real, n_immersion::Real)
    # Calculate interface-specific coefficients
    Tp_mc, Ts_mc = calculate_interface_fresnel(kr2, λ, n_medium, n_coverslip)
    Tp_ci, Ts_ci = calculate_interface_fresnel(kr2, λ, n_coverslip, n_immersion)

    # Calculate combined transmission coefficients
    Tp = Tp_mc * Tp_ci
    Ts = Ts_mc * Ts_ci

    return Tp, Ts
end

"""
    calculate_axial_phase(z::Real, z_stage::Real, kz_medium::Complex, 
                         kz_coverslip::Complex, kz_immersion::Complex,
                         coverslip_thickness::Real=0.17)

Calculate total axial phase from defocus.

# Arguments
- `z`: Emitter position relative to the coverslip (depth above the coverslip in μm)
- `z_stage`: Distance the sample stage was moved away from the nominal focal plane at the coverslip (μm)
- `kz_medium`: z-component of wave vector in sample medium
- `kz_coverslip`: z-component of wave vector in coverslip
- `kz_immersion`: z-component of wave vector in immersion medium
- `coverslip_thickness`: Thickness of coverslip in mm (default: 0.17mm)

# Returns
- Defocus phase

"""
function calculate_axial_phase(z::Real, z_stage::Real, kz_medium::Complex,
    kz_coverslip::Complex, kz_immersion::Complex,
    coverslip_thickness::Real=0.17)
    # Emitter position relative to coverslip interface

    return kz_medium * z - kz_immersion * z_stage
end

"""
    calculate_dipole_field_components(ϕ::Complex, sinθ::Complex, cosθ::Complex,
                                    dipole::DipoleVector, Tp::Complex, Ts::Complex)

Calculate vectorial field components for a dipole orientation.

# Arguments
- `ϕ`: Azimuthal angle in pupil plane
- `sinθ`: Sine of polar angle
- `cosθ`: Cosine of polar angle
- `dipole`: Dipole orientation vector
- `Tp`: p-polarization transmission coefficient
- `Ts`: s-polarization transmission coefficient

# Returns
- Tuple (Ex, Ey) of complex field components
"""
function calculate_dipole_field_components(ϕ::Complex, sinθ::Complex, cosθ::Complex,
    dipole::DipoleVector, Tp::Complex, Ts::Complex)
    # Compute the vectorial field components based on the dipole orientation
    # For a dipole, the field has both transverse and longitudinal components
    Eθ = (Tp * (dipole.px * cos(ϕ) + dipole.py * sin(ϕ)) * cosθ -
          Tp * dipole.pz * sinθ)
    Eϕ = Ts * (-dipole.px * sin(ϕ) + dipole.py * cos(ϕ))

    # Convert to Cartesian coordinates
    Ex = Eθ * cos(ϕ) - Eϕ * sin(ϕ)
    Ey = Eθ * sin(ϕ) + Eϕ * cos(ϕ)

    return Ex, Ey
end

"""
    calculate_apodization(kr2::Real, λ::Real, 
                        n_medium::Real, n_immersion::Real)

Calculate apodization factor for energy conservation.

# Arguments
- `kr2`: Squared lateral spatial frequency
- `λ`: Wavelength in microns
- `n_medium`: Refractive index of the sample medium
- `n_immersion`: Refractive index of the immersion medium

# Returns
- Apodization factor for energy conservation
"""
function calculate_apodization(kr2::Real, λ::Real,
    n_medium::Real, n_immersion::Real)
    # Calculate cosines of angles
    cosθ_medium = sqrt(complex(1 - kr2 * λ^2 / n_medium^2))
    cosθ_immersion = sqrt(complex(1 - kr2 * λ^2 / n_immersion^2))

    # Calculate apodization factor
    # This accounts for the change in solid angle with refraction
    apod = sqrt(cosθ_immersion) / cosθ_medium

    return apod
end
