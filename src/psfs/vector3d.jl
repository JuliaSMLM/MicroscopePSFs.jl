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
- Vector3DPSF instance
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
    
    # Store the base pupil and/or convert Zernike coeffs to pupil if provided
    stored_base_pupil = base_pupil
    stored_zernike = base_zernike
    
    if isnothing(base_pupil) && !isnothing(base_zernike)
        # Create base pupil from Zernike coefficients
        base_pupil = PupilFunction(nₐ, λ, n_medium, base_zernike; grid_size=grid_size)
    end
    
    # Create vector pupil
    vpupil = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, grid_size)
    
    # Fill with vector components using base aberration
    if !isnothing(base_pupil)
        fill_vector_pupil!(vpupil, dipole, focal_z, base_pupil)
    else
        fill_vector_pupil!(vpupil, dipole, focal_z)
    end
    
    return Vector3DPSF{T}(
        T(nₐ), T(λ), T(n_medium), T(n_coverslip),
        T(n_immersion), dipole, T(focal_z), vpupil,
        stored_base_pupil, stored_zernike
    )
end

"""
    amplitude(psf::Vector3DPSF, x::Real, y::Real, z::Real)

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
function amplitude(psf::Vector3DPSF, x::Real, y::Real, z::Real)
    return amplitude(psf.pupil, x, y, z)
end

"""
    (psf::Vector3DPSF)(x::Real, y::Real, z::Real)

Compute PSF intensity at given position by summing squared field amplitudes.

# Arguments
- `x, y`: Lateral position in microns relative to PSF center
- `z`: Axial position in microns relative to focal plane

# Returns
- Total intensity |Ex|² + |Ey|²

# Notes
- For randomly oriented dipoles, evaluate three orthogonal orientations and average
"""
function (psf::Vector3DPSF)(x::Real, y::Real, z::Real)
    E = amplitude(psf, x, y, z)
    return abs2(E[1]) + abs2(E[2])
end

"""
    update_pupil!(psf::Vector3DPSF) -> Vector3DPSF

Update the vector pupil function based on stored Zernike coefficients and/or base pupil.
This is useful after modifying aberrations to regenerate the pupil fields.

# Arguments
- `psf`: Vector3DPSF to update

# Returns
- Updated Vector3DPSF

# Notes
- Requires either stored Zernike coefficients or a base pupil
- Returns the updated PSF for method chaining
"""
function update_pupil!(psf::Vector3DPSF)
    # Check if we have necessary components to update
    if isnothing(psf.base_pupil) && isnothing(psf.zernike_coeffs)
        throw(ArgumentError("Cannot update pupil: no base pupil or Zernike coefficients stored"))
    end
    
    # If we have Zernike coefficients but no base pupil, create it
    updated_base_pupil = psf.base_pupil
    if isnothing(updated_base_pupil) && !isnothing(psf.zernike_coeffs)
        updated_base_pupil = PupilFunction(
            psf.nₐ, psf.λ, psf.n_medium, 
            psf.zernike_coeffs; 
            grid_size=size(psf.pupil.Ex.field, 1)
        )
    end
    
    # Fill the vector pupil with updated components
    if !isnothing(updated_base_pupil)
        fill_vector_pupil!(psf.pupil, psf.dipole, psf.focal_z, updated_base_pupil)
    else
        fill_vector_pupil!(psf.pupil, psf.dipole, psf.focal_z)
    end
    
    return psf
end

"""
    integrate_pixels(psf::Vector3DPSF,
                    camera::AbstractCamera,
                    emitter::AbstractEmitter;
                    sampling::Integer=2)

Integrate Vector3DPSF over camera pixels.

# Arguments
- `psf`: Vector3DPSF instance with fixed dipole orientation
- `camera`: Camera geometry
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density for integration accuracy

# Returns
- Array of integrated PSF intensities with dimensions [ny, nx]
- Values represent actual photon counts based on emitter.photons

# Notes
- Dipole orientation comes from the PSF itself, not the emitter
- For varying dipole orientations, create multiple PSFs
"""
function integrate_pixels(
    psf::Vector3DPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    # Check if emitter has required z-coordinate
    if !hasfield(typeof(emitter), :z)
        throw(ArgumentError("Vector3DPSF requires an emitter with a z-coordinate"))
    end
    
    # Use the generic integration method
    result = _integrate_pixels_generic(
        psf, camera, emitter,
        (p, x, y) -> p(x, y, emitter.z),
        Float64; sampling=sampling
    )
    
    # Multiply by photon count to preserve physical meaning
    return result .* emitter.photons
end

"""
    integrate_pixels_amplitude(psf::Vector3DPSF,
                              camera::AbstractCamera,
                              emitter::AbstractEmitter;
                              sampling::Integer=2)

Integrate Vector3DPSF complex amplitude over camera pixels.

# Arguments
- `psf`: Vector3DPSF instance with fixed dipole orientation
- `camera`: Camera geometry
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density for integration accuracy

# Returns
- Array of integrated complex field components with dimensions [ny, nx, 2]
- First two dimensions are spatial, third dimension holds [Ex, Ey]

# Notes
- For coherent calculations in vectorial microscopy
- Preserves relative phase information between field components
"""
function integrate_pixels_amplitude(
    psf::Vector3DPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    # Check if emitter has required z-coordinate
    if !hasfield(typeof(emitter), :z)
        throw(ArgumentError("Vector3DPSF requires an emitter with a z-coordinate"))
    end
    
    # Integration for complex field components
    Ex = _integrate_pixels_generic(
        psf, camera, emitter,
        (p, x, y) -> amplitude(p, x, y, emitter.z)[1],
        Complex{Float64}; sampling=sampling
    )
    
    Ey = _integrate_pixels_generic(
        psf, camera, emitter,
        (p, x, y) -> amplitude(p, x, y, emitter.z)[2],
        Complex{Float64}; sampling=sampling
    )
    
    # Combine into a 3D array: [ny, nx, 2] for Ex, Ey components
    result = Array{Complex{Float64}, 3}(undef, size(Ex,1), size(Ex,2), 2)
    result[:,:,1] = Ex
    result[:,:,2] = Ey
    
    return result
end

# Display method
function Base.show(io::IO, psf::Vector3DPSF)
    print(io, "Vector3DPSF(NA=$(psf.nₐ), λ=$(psf.λ)μm, n_medium=$(psf.n_medium))")
    has_bp = !isnothing(psf.base_pupil)
    has_zc = !isnothing(psf.zernike_coeffs)
    
    if has_bp || has_zc
        print(io, " with ")
        if has_bp
            print(io, "base pupil")
            has_zc && print(io, " and ")
        end
        if has_zc
            print(io, "$(length(psf.zernike_coeffs)) Zernike terms")
        end
    end
end