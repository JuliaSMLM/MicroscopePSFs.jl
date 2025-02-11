# src/psfs/vector3d.jl

"""
    calculate_fresnel_coefficients(kr2::Real, λ::Real, 
                                 n_medium::Real, n_immersion::Real)

Compute Fresnel transmission coefficients including super-critical angle fluorescence.

# Arguments
- `kr2`: Radial spatial frequency squared
- `λ`: Wavelength in microns
- `n_medium`: Sample medium refractive index
- `n_immersion`: Immersion medium refractive index

# Returns
- Tuple of (Tp, Ts, sinθ₁, cosθ₁, apod):
  * Tp: p-polarization transmission coefficient
  * Ts: s-polarization transmission coefficient
  * sinθ₁: Sine of angle in medium
  * cosθ₁: Cosine of angle in medium
  * apod: Apodization factor for angular spectrum
"""
function calculate_fresnel_coefficients(kr2::Real, λ::Real, 
                                     n_medium::Real, n_immersion::Real)
    # Calculate angles using complex sqrt for automatic handling of evanescent waves
    sinθ₁ = complex(sqrt(kr2)*λ/n_medium)  # Explicitly make complex
    cosθ₁ = sqrt(complex(1 - kr2*λ^2/n_medium^2))
    cosθᵢ = sqrt(complex(1 - kr2*λ^2/n_immersion^2))
    
    # Calculate Fresnel coefficients
    Tp = 2.0*n_medium*cosθ₁/(n_medium*cosθᵢ + n_immersion*cosθ₁)
    Ts = 2.0*n_medium*cosθ₁/(n_medium*cosθ₁ + n_immersion*cosθᵢ)
    
    # Apodization factor including SAF contribution
    apod = sqrt(cosθᵢ)/cosθ₁
    
    return Tp, Ts, sinθ₁, cosθ₁, apod
end

"""
    calculate_wave_vectors(kr2::Real, λ::Real, n_medium::Real, n_immersion::Real)

Calculate wave vectors for both media including evanescent waves.

# Arguments
- `kr2`: Radial spatial frequency squared
- `λ`: Wavelength in microns
- `n_medium`: Sample medium refractive index
- `n_immersion`: Immersion medium refractive index

# Returns
- Tuple of (kz_medium, kz_immersion): z-components of wave vectors in each medium
"""
function calculate_wave_vectors(kr2::Real, λ::Real, n_medium::Real, n_immersion::Real)
    k₀ = 2π/λ
    
    # Wave vector z-components in both media (can be complex)
    kz_medium = k₀ * sqrt(complex(n_medium^2 - kr2*λ^2))
    kz_immersion = k₀ * sqrt(complex(n_immersion^2 - kr2*λ^2))
    
    return kz_medium, kz_immersion
end

"""
    calculate_field_components(ϕ::Real, Tp::Complex, Ts::Complex, 
                             sinθ₁::Complex, cosθ₁::Complex, 
                             dipole::DipoleVector)

Compute vectorial field components including polarization effects.

# Arguments
- `ϕ`: Azimuthal angle in pupil plane
- `Tp`: p-polarization transmission coefficient
- `Ts`: s-polarization transmission coefficient
- `sinθ₁`: Sine of angle in medium
- `cosθ₁`: Cosine of angle in medium
- `dipole`: Dipole orientation vector

# Returns
- Tuple (Ex, Ey) of complex field components
"""
function calculate_field_components(ϕ::Complex, Tp::Complex, Ts::Complex,
                                 sinθ₁::Complex, cosθ₁::Complex,
                                 dipole::DipoleVector)
    # Compute the vectorial field components based on the dipole orientation
    Eθ = (Tp * (dipole.px*cos(ϕ) + dipole.py*sin(ϕ)) * cosθ₁ - 
          Tp * dipole.pz * sinθ₁)
    Eϕ = Ts * (-dipole.px*sin(ϕ) + dipole.py*cos(ϕ))
    
    # Convert to Cartesian coordinates
    Ex = Eθ * cos(ϕ) - Eϕ * sin(ϕ)
    Ey = Eθ * sin(ϕ) + Eϕ * cos(ϕ)
    
    return Ex, Ey
end

"""
    Vector3DPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
                base_pupil::Union{Nothing, PupilFunction}=nothing,
                base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
                n_medium::Real=1.33,
                n_coverslip::Real=1.52,
                n_immersion::Real=1.52,
                focal_z::Real=0.0,
                grid_size::Integer=128)

Create a vector PSF using either a pupil-based or Zernike-based approach.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `dipole`: Dipole orientation vector

# Keyword Arguments
- `base_pupil`: Optional base aberration pupil function
- `base_zernike`: Optional Zernike coefficients for base aberration
- `n_medium`: Sample medium refractive index (default: 1.33)
- `n_coverslip`: Cover slip refractive index (default: 1.52)
- `n_immersion`: Immersion medium refractive index (default: 1.52)
- `focal_z`: Focal plane position in microns (default: 0.0)
- `grid_size`: Size of pupil grid (default: 128)

# Returns
- Vector3DPupilPSF instance
"""
function Vector3DPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
    base_pupil::Union{Nothing, PupilFunction}=nothing,
    base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
    n_medium::Real=1.33,
    n_coverslip::Real=1.52,
    n_immersion::Real=1.52,
    focal_z::Real=0.0,
    grid_size::Integer=128)
    
    T = promote_type(typeof(nₐ), typeof(λ), typeof(n_medium), typeof(focal_z))
    
    # Convert Zernike coeffs to pupil if provided
    base_pupil = if !isnothing(base_zernike)
        PupilFunction(nₐ, λ, n_medium, base_zernike; grid_size=grid_size)
    else
        base_pupil
    end
    
    # Create vector pupil
    vpupil = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, grid_size)
    
    # Fill with vector components using base aberration
    if !isnothing(base_pupil)
        fill_vector_pupil!(vpupil, dipole, focal_z, base_pupil)
    else
        fill_vector_pupil!(vpupil, dipole, focal_z)
    end
    
    return Vector3DPupilPSF{T}(
        T(nₐ), T(λ), T(n_medium), T(n_coverslip),
        T(n_immersion), dipole, T(focal_z), vpupil
    )
end

"""
    amplitude(psf::Vector3DPupilPSF, x::Real, y::Real, z::Real)

Compute complex vector amplitude at given position.

# Arguments
- `psf`: Vector PSF instance
- `x, y`: Lateral position in microns relative to PSF center
- `z`: Axial position in microns relative to focal plane

# Returns
- Vector [Ex, Ey] of complex field amplitudes

# Notes
- z coordinate is relative to current focal plane position (psf.focal_z)
- Includes both UAF and SAF contributions automatically
"""
function amplitude(psf::Vector3DPupilPSF, x::Real, y::Real, z::Real)
    return amplitude(psf.pupil, x, y, z)
end

"""
    (psf::Vector3DPupilPSF)(x::Real, y::Real, z::Real)

Compute PSF intensity at given position by summing squared field amplitudes.

# Arguments
- `x, y`: Lateral position in microns relative to PSF center
- `z`: Axial position in microns relative to focal plane

# Returns
- Total intensity |Ex|² + |Ey|²

# Notes
- For randomly oriented dipoles, evaluate three orthogonal orientations and average
"""
function (psf::Vector3DPupilPSF)(x::Real, y::Real, z::Real)
    E = amplitude(psf, x, y, z)
    return abs2(E[1]) + abs2(E[2])
end