# src/psfs/vector3d.jl

"""
    Vector3DPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
                base_pupil::Union{Nothing, PupilFunction}=nothing,
                base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
                n_medium::Real=1.33,
                n_coverslip::Real=1.52,
                n_immersion::Real=1.52,
                z_stage::Real=0.0,
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
- `z_stage`: Distance the sample stage was moved away from the nominal focal plane at the coverslip (μm) (default: 0.0)
- `grid_size`: Size of pupil grid (default: 128)

# Returns
- Vector3DPSF instance
"""
function Vector3DPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
    base_pupil::Union{Nothing,PupilFunction}=nothing,
    base_zernike::Union{Nothing,ZernikeCoefficients}=nothing,
    n_medium::Real=1.33,
    n_coverslip::Real=1.52,
    n_immersion::Real=1.52,
    z_stage::Real=0.0,
    grid_size::Integer=128)

    T = promote_type(typeof(nₐ), typeof(λ), typeof(n_medium), typeof(z_stage))

    # Store the base pupil and/or convert Zernike coeffs to pupil if provided
    stored_base_pupil = base_pupil
    stored_zernike = base_zernike

    if isnothing(base_pupil) && !isnothing(base_zernike)
        # Create base pupil from Zernike coefficients
        base_pupil = PupilFunction(nₐ, λ, n_medium, base_zernike; grid_size=grid_size)
    end

    # Create vector pupil function with all position-independent factors
    vector_pupils = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, grid_size)

    # Fill with vector components using base aberration
    if !isnothing(base_pupil)
        fill_vector_pupils!(vector_pupils, dipole, base_pupil)
    else
        fill_vector_pupils!(vector_pupils, dipole)
    end

    return Vector3DPSF{T}(
        T(nₐ), T(λ), T(n_medium), T(n_coverslip),
        T(n_immersion), dipole, T(z_stage), vector_pupils,
        stored_base_pupil, stored_zernike
    )
end

"""
    amplitude(psf::Vector3DPSF, x::Real, y::Real, z::Real)

Compute complex vector amplitude at given position.

# Arguments
- `psf`: Vector PSF instance
- `x, y`: Lateral position in microns relative to PSF center
- `z`: Axial position in microns representing depth above the coverslip

# Returns
- Vector [Ex, Ey] of complex field amplitudes

# Notes
- z coordinate represents the depth above the coverslip
- z_stage in the PSF indicates the distance the stage was moved away from the nominal focal plane
- Includes both UAF and SAF contributions automatically
"""
function amplitude(psf::Vector3DPSF{T}, x::Real, y::Real, z::Real) where {T}

    # Initialize result vector [Ex, Ey]
    result = zeros(Complex{promote_type(T, typeof(x), typeof(y), typeof(z))}, 2)

    # Pupil parameters
    vector_pupils = psf.vector_pupils
    grid_size = size(vector_pupils.Ex.field, 1)
    kmax_val = psf.nₐ / psf.λ
    kpixel = 2 * kmax_val / (grid_size - 1)
    center = (grid_size + 1) / 2

    # Integrate over pupil
    for i in 1:grid_size, j in 1:grid_size
        # Get k-space coordinates
        kx = (i - center) * kpixel
        ky = (j - center) * kpixel
        kr2 = kx^2 + ky^2

        # Skip points outside the pupil
        kr2 > kmax_val^2 && continue

        # Get wave vectors in each medium
        kz_medium, kz_coverslip, kz_immersion = calculate_wave_vectors(
            kr2, psf.λ, psf.n_medium, psf.n_coverslip, psf.n_immersion)

        # Calculate total phase with position dependence
        lateral_phase = kx * x + ky * y
        axial_phase = calculate_axial_phase(z, psf.z_stage, kz_medium, kz_coverslip, kz_immersion)
        total_phase = 2π * (lateral_phase + axial_phase)

        # Apply phase to pre-calculated pupil field
        phase_factor = exp(im * total_phase)

        # Add contribution to result
        result[1] += vector_pupils.Ex.field[j, i] * phase_factor * kpixel^2
        result[2] += vector_pupils.Ey.field[j, i] * phase_factor * kpixel^2
    end

    return result
end

"""
    (psf::Vector3DPSF)(x::Real, y::Real, z::Real)

Compute PSF intensity at given position by summing squared field amplitudes.

# Arguments
- `x, y`: Lateral position in microns relative to PSF center
- `z`: Axial position in microns representing depth above the coverslip

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
    update_pupils!(psf::Vector3DPSF) -> Vector3DPSF

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
function update_pupils!(psf::Vector3DPSF)
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
            grid_size=size(psf.vector_pupils.Ex.field, 1)
        )
    end

    # Fill the vector pupil with updated components
    if !isnothing(updated_base_pupil)
        fill_vector_pupils!(psf.vector_pupils, psf.dipole, updated_base_pupil)
    else
        fill_vector_pupils!(psf.vector_pupils, psf.dipole)
    end

    return psf
end


"""
    integrate_pixels(
        psf::Vector3DPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2
    )

Integrate Vector3DPSF over camera pixels.

# Arguments
- `psf`: Vector3DPSF instance with fixed dipole orientation
- `camera`: Camera geometry
- `emitter`: Emitter with position information (emitter.z represents depth above the coverslip)
- `sampling`: Subpixel sampling density for integration accuracy

# Returns
- Array of integrated PSF intensities with dimensions [ny, nx]
- Values represent actual photon counts based on emitter.photons

# Notes
- Dipole orientation comes from the PSF itself, not the emitter
- For varying dipole orientations, create multiple PSFs
- The emitter's z position must represent depth above the coverslip
"""
function integrate_pixels(
    psf::Vector3DPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    # Check if emitter has required z-coordinate and provide domain-specific error message
    if !has_z_coordinate(emitter)
        throw(ArgumentError("Vector3DPSF requires an emitter with a z-coordinate (depth above the coverslip)"))
    end

    # Use the standard integration function - it will automatically handle the z-coordinate
    return invoke(
        integrate_pixels,
        Tuple{AbstractPSF,AbstractCamera,AbstractEmitter},
        psf, camera, emitter;
        sampling=sampling
    )
end

"""
    integrate_pixels_amplitude(
        psf::Vector3DPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2
    )

Integrate Vector3DPSF complex amplitude over camera pixels.

# Arguments
- `psf`: Vector3DPSF instance with fixed dipole orientation
- `camera`: Camera geometry
- `emitter`: Emitter with position information (emitter.z represents depth above the coverslip)
- `sampling`: Subpixel sampling density for integration accuracy

# Returns
- Array of integrated complex field components with dimensions [ny, nx, 2]
- First two dimensions are spatial, third dimension holds [Ex, Ey]

# Notes
- For coherent calculations in vectorial microscopy
- Preserves relative phase information between field components
- The emitter's z position must represent depth above the coverslip
"""
function integrate_pixels_amplitude(
    psf::Vector3DPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    # Check if emitter has required z-coordinate and provide domain-specific error message
    if !has_z_coordinate(emitter)
        throw(ArgumentError("Vector3DPSF requires an emitter with a z-coordinate (depth above the coverslip)"))
    end

    # Use the standard integration function - it will automatically handle the z-coordinate
    return invoke(
        integrate_pixels_amplitude,
        Tuple{AbstractPSF,AbstractCamera,AbstractEmitter},
        psf, camera, emitter;
        sampling=sampling
    )
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
