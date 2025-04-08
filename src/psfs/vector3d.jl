# src/psfs/vector3d.jl

"""
    VectorPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
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
- VectorPSF instance
"""
function VectorPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
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

    # Create a single vector pupil function
    vector_pupil = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, grid_size)

    # Fill with vector components using base aberration
    if !isnothing(base_pupil)
        fill_vector_pupils!(vector_pupil, dipole, base_pupil; normalize=true)
    else
        fill_vector_pupils!(vector_pupil, dipole; normalize=true)
    end

    # Create vector with single pupil
    vector_pupils = [vector_pupil]

    return VectorPSF{T}(
        T(nₐ), T(λ), T(n_medium), T(n_coverslip),
        T(n_immersion), dipole, T(z_stage), vector_pupils,
        stored_base_pupil, stored_zernike
    )
end

"""
    VectorPSF(nₐ::Real, λ::Real;
                base_pupil::Union{Nothing, PupilFunction}=nothing,
                base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
                n_medium::Real=1.33,
                n_coverslip::Real=1.52,
                n_immersion::Real=1.52,
                z_stage::Real=0.0,
                grid_size::Integer=128)

Create a vector PSF with a rotating dipole, modeled as an incoherent sum of 
x, y, and z dipole orientations.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns

# Keyword Arguments
- `base_pupil`: Optional base aberration pupil function
- `base_zernike`: Optional Zernike coefficients for base aberration
- `n_medium`: Sample medium refractive index (default: 1.33)
- `n_coverslip`: Cover slip refractive index (default: 1.52)
- `n_immersion`: Immersion medium refractive index (default: 1.52)
- `z_stage`: Distance the sample stage was moved away from the nominal focal plane at the coverslip (μm) (default: 0.0)
- `grid_size`: Size of pupil grid (default: 128)

# Returns
- VectorPSF instance representing a rotating dipole

# Notes
- Intensity is calculated as incoherent average of three orthogonal dipole orientations
- The amplitude() function is not meaningful for rotating dipoles and will throw an error
"""
function VectorPSF(nₐ::Real, λ::Real;
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

    # Create dipole orientations for x, y, z
    dipoles = [
        DipoleVector(1.0, 0.0, 0.0),  # x-dipole
        DipoleVector(0.0, 1.0, 0.0),  # y-dipole
        DipoleVector(0.0, 0.0, 1.0)   # z-dipole
    ]
    
    # Use a placeholder dipole for the main struct field
    # This isn't used for calculations with multiple pupils, just for metadata
    placeholder_dipole = DipoleVector(0.0, 0.0, 0.0)
    
    # Create vector pupils for each dipole orientation
    vector_pupils = Vector{VectorPupilFunction{T}}(undef, 3)
    
    vector_pupils = Vector{VectorPupilFunction{T}}(undef, 3)
    
    for i in 1:3
        pupil = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, grid_size)
        if !isnothing(base_pupil)
            fill_vector_pupils!(pupil, dipoles[i], base_pupil, normalize=false)
        else
            fill_vector_pupils!(pupil, dipoles[i], normalize=false)
        end
        vector_pupils[i] = pupil
    end
    
    # Calculate total energy across all pupils
    total_energy = 0.0
    for pupil in vector_pupils
        kpix² = kpixelsize(pupil)^2
        total_energy += (sum(abs2, pupil.Ex.field) + sum(abs2, pupil.Ey.field)) * kpix²
    end
    
    # Apply joint normalization
    scale = 1/sqrt(total_energy)
    for pupil in vector_pupils
        pupil.Ex.field .*= scale
        pupil.Ey.field .*= scale
    end

    return VectorPSF{T}(
        T(nₐ), T(λ), T(n_medium), T(n_coverslip),
        T(n_immersion), placeholder_dipole, T(z_stage), vector_pupils,
        stored_base_pupil, stored_zernike
    )
end

"""
    _calculate_field_amplitude(
        pupil::VectorPupilFunction,
        nₐ::Real, λ::Real, n_medium::Real, n_coverslip::Real, n_immersion::Real, z_stage::Real,
        x::Real, y::Real, z::Real)

Internal helper function to compute complex vector amplitude at given position.

# Arguments
- `pupil`: Vector pupil function containing Ex and Ey fields
- `nₐ`, `λ`, `n_medium`, `n_coverslip`, `n_immersion`: Optical parameters 
- `z_stage`: Distance the sample stage was moved away from the nominal focal plane
- `x, y, z`: Position in microns

# Returns
- Vector [Ex, Ey] of complex field amplitudes

# Notes
- This is a helper function used by both amplitude() and the PSF evaluation function
- Implements the core field calculation for a single pupil
"""
function _calculate_field_amplitude(
    pupil::VectorPupilFunction,
    nₐ::Real, λ::Real, n_medium::Real, n_coverslip::Real, n_immersion::Real, z_stage::Real,
    x::Real, y::Real, z::Real)
    
    # Initialize result vector [Ex, Ey]
    # Fixed: Don't use the type constructor in this way
    T = promote_type(real(eltype(pupil.Ex.field)), typeof(x), typeof(y), typeof(z))
    result = zeros(Complex{T}, 2)

    # Pupil parameters
    grid_size = size(pupil.Ex.field, 1)
    kmax_val = nₐ / λ
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
            kr2, λ, n_medium, n_coverslip, n_immersion)

        # Calculate total phase with position dependence
        lateral_phase = kx * x + ky * y
        axial_phase = calculate_axial_phase(z, z_stage, kz_medium, kz_coverslip, kz_immersion)
        total_phase = 2π * (lateral_phase + axial_phase)

        # Apply phase to pre-calculated pupil field
        phase_factor = exp(im * total_phase)

        # Add contribution to result
        result[1] += pupil.Ex.field[j, i] * phase_factor * kpixel^2
        result[2] += pupil.Ey.field[j, i] * phase_factor * kpixel^2
    end

    return result
end

"""
    amplitude(psf::VectorPSF, x::Real, y::Real, z::Real)

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
- Not meaningful for rotating dipoles (throws an error if PSF has multiple pupils)
"""
function amplitude(psf::VectorPSF{T}, x::Real, y::Real, z::Real) where {T}
    if length(psf.vector_pupils) > 1
        error("amplitude() not meaningful for multiple dipole orientations; use intensity evaluation directly")
    end
    
    # Use the first (and only) pupil
    return _calculate_field_amplitude(
        psf.vector_pupils[1],
        psf.nₐ, psf.λ, psf.n_medium, psf.n_coverslip, psf.n_immersion, psf.z_stage,
        x, y, z
    )
end

"""
    (psf::VectorPSF)(x::Real, y::Real, z::Real)

Compute PSF intensity at given position by summing squared field amplitudes.

# Arguments
- `x, y`: Lateral position in microns relative to PSF center
- `z`: Axial position in microns representing depth above the coverslip

# Returns
- Total intensity |Ex|² + |Ey|²

# Notes
- For single dipole: evaluates intensity from field amplitude
- For rotating dipole: calculates incoherent sum over three orthogonal dipole orientations
"""
function (psf::VectorPSF)(x::Real, y::Real, z::Real)
    if length(psf.vector_pupils) == 1
        # Single dipole case
        E = _calculate_field_amplitude(
            psf.vector_pupils[1],
            psf.nₐ, psf.λ, psf.n_medium, psf.n_coverslip, psf.n_immersion, psf.z_stage,
            x, y, z
        )
        return abs2(E[1]) + abs2(E[2])
    else
        # Multiple dipoles case - incoherent sum
        intensity = 0.0
        
        for pupil in psf.vector_pupils
            E = _calculate_field_amplitude(
                pupil,
                psf.nₐ, psf.λ, psf.n_medium, psf.n_coverslip, psf.n_immersion, psf.z_stage,
                x, y, z
            )
            intensity += abs2(E[1]) + abs2(E[2])
        end
        
        # Return average intensity
        return intensity 
    end
end

"""
    update_pupils!(psf::VectorPSF) -> VectorPSF

Update the vector pupil functions based on stored Zernike coefficients and/or base pupil.
This is useful after modifying aberrations to regenerate the pupil fields.

# Arguments
- `psf`: VectorPSF to update

# Returns
- Updated VectorPSF

# Notes
- Requires either stored Zernike coefficients or a base pupil
- For rotating dipoles (multiple pupils), updates each pupil with the appropriate dipole orientation
- Returns the updated PSF for method chaining
"""
function update_pupils!(psf::VectorPSF)
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
            grid_size=size(psf.vector_pupils[1].Ex.field, 1)
        )
    end

    # For single dipole case
    if length(psf.vector_pupils) == 1
        # Fill the vector pupil with updated components
        if !isnothing(updated_base_pupil)
            fill_vector_pupils!(psf.vector_pupils[1], psf.dipole, updated_base_pupil; normalize=true)
        else
            fill_vector_pupils!(psf.vector_pupils[1], psf.dipole; normalize=true)
        end
    else
        # For rotating dipole case, update each pupil with its corresponding dipole orientation
        dipoles = [
            DipoleVector(1.0, 0.0, 0.0),  # x-dipole
            DipoleVector(0.0, 1.0, 0.0),  # y-dipole
            DipoleVector(0.0, 0.0, 1.0)   # z-dipole
        ]
        
        for i in 1:length(psf.vector_pupils)
            if !isnothing(updated_base_pupil)
                fill_vector_pupils!(psf.vector_pupils[i], dipoles[i], updated_base_pupil)
            else
                fill_vector_pupils!(psf.vector_pupils[i], dipoles[i])
            end
        end
    end

    return psf
end

"""
    integrate_pixels_amplitude!(
        result::AbstractArray{Complex{T},3},
        psf::VectorPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2,
        threaded::Bool=true
    ) where T <: Real

Special version of amplitude integration for VectorPSF that preserves polarization components.
Uses the flexible core integration function that handles vector returns.

# Arguments
- `result`: Pre-allocated 3D complex array where results will be stored, dimensions [y, x, pol]
- `psf`: VectorPSF instance
- `camera`: Camera geometry defining pixel edges
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- The `result` array, filled with integrated complex amplitudes for each polarization
"""
function integrate_pixels_amplitude!(
    result::AbstractArray{Complex{T},3},
    psf::VectorPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2,
    threaded::Bool=true
) where T <: Real
    # Verify dimensions
    ny, nx, npol = size(result)
    npol == 2 || throw(DimensionMismatch("Third dimension of result array must be 2 for Ex and Ey"))
    
    # Use the core integration function directly - it now handles vector returns
    _integrate_pixels_generic!(
        result,
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitter,
        amplitude,  # The amplitude function returns [Ex, Ey]
        sampling=sampling,
        threaded=threaded
    )
    
    return result
end

"""
    integrate_pixels_amplitude(
        psf::VectorPSF,
        camera::AbstractCamera,
        emitter::AbstractEmitter;
        sampling::Integer=2,
        threaded::Bool=true
    )

Specialized version of amplitude integration for VectorPSF that preserves polarization components.

# Arguments
- `psf`: VectorPSF instance
- `camera`: Camera geometry defining pixel edges
- `emitter`: Emitter with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- 3D array of integrated complex amplitudes with dimensions [y, x, pol]
  where pol index 1 = Ex and pol index 2 = Ey
"""
function integrate_pixels_amplitude(
    psf::VectorPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2,
    threaded::Bool=true
)
    if length(psf.vector_pupils) > 1
        error("integrate_pixels_amplitude() not meaningful for multiple dipole orientations")
    end
    
    T = Complex{typeof(emitter.photons)}
    
    # Use the generic function with additional dimension specified
    return _integrate_pixels_generic(
        psf,
        camera.pixel_edges_x,
        camera.pixel_edges_y,
        emitter,
        amplitude,
        T,
        (2,);  # Additional dimension for the polarization components
        sampling=sampling,
        threaded=threaded
    )
end

"""
    integrate_pixels_amplitude(
        psf::VectorPSF,
        camera::AbstractCamera,
        emitters::Vector{<:AbstractEmitter};
        sampling::Integer=2,
        threaded::Bool=true
    )

Integrate PSF complex amplitude over camera pixels for multiple emitters.
Special version for VectorPSF that preserves polarization components.

# Arguments
- `psf`: VectorPSF instance
- `camera`: Camera geometry defining pixel edges
- `emitters`: Vector of emitters with position information
- `sampling`: Subpixel sampling density (default: 2)
- `threaded`: Whether to use multi-threading for integration (default: true)

# Returns
- 3D array of integrated complex amplitudes with dimensions [y, x, pol]
  where pol index 1 = Ex and pol index 2 = Ey

# Notes
- Results are the coherent sum of individual emitter field contributions
"""
function integrate_pixels_amplitude(
    psf::VectorPSF,
    camera::AbstractCamera,
    emitters::Vector{<:AbstractEmitter};
    sampling::Integer=2,
    threaded::Bool=true
)
    if length(psf.vector_pupils) > 1
        error("integrate_pixels_amplitude() not meaningful for multiple dipole orientations")
    end
    
    isempty(emitters) && return zeros(Complex{Float64}, length(camera.pixel_edges_y)-1, length(camera.pixel_edges_x)-1, 2)
    
    # Determine result type from first emitter
    T = Complex{typeof(emitters[1].photons)}
    
    # Allocate result array
    ny = length(camera.pixel_edges_y) - 1
    nx = length(camera.pixel_edges_x) - 1
    result = zeros(T, ny, nx, 2)  # 3D array for polarization components
    
    # Temporary buffer for individual emitter contribution
    buffer = similar(result)
    
    # Process each emitter
    for emitter in emitters
        fill!(buffer, zero(T))
        integrate_pixels_amplitude!(buffer, psf, camera, emitter; sampling=sampling, threaded=threaded)
        result .+= buffer
    end
    
    return result
end