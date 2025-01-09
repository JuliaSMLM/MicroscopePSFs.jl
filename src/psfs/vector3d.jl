"""
    calculate_pupil_field(ϕ, Tp, Ts, sinθ₁, cosθ₁, dipole::DipoleVector)

Computes electric field components in the pupil plane for a dipole emitter.

# Physics
- Uses angular spectrum representation for high-NA imaging
- Includes both p-(radial) and s-(azimuthal) polarization components
- Accounts for Fresnel transmission coefficients (Tp, Ts)
- Optional phase retardation δ for birefringent media

# Parameters
- ϕ: Azimuthal angle in pupil plane
- Tp: p-polarization Fresnel transmission coefficient
- Ts: s-polarization Fresnel transmission coefficient
- sinθ₁, cosθ₁: Sine/cosine of polar angle in medium 1
- δ: Phase retardation (optional)
- dvec: Dipole moment orientation vector [default: [1,1,1]]

# Returns
- Ex, Ey: Electric field components in x,y
- field_matrix: Matrix of field components
"""
function calculate_pupil_field(ϕ, Tp, Ts, sinθ₁, cosθ₁, dipole::DipoleVector)
    # Normalize dipole vector
    dvec = [dipole.px, dipole.py, dipole.pz] ./ norm([dipole.px, dipole.py, dipole.pz])
    
    # Compute the vectorial field components based on the dipole orientation
    Eθ = (Tp .* (dvec[1]*cosϕ + dvec[2]*sinϕ) * cosθ₁ - Tp .* dvec[3] * sinθ₁)
    Eϕ = Ts .* (-dvec[1]*sinϕ + dvec[2]*cosϕ)
    
    # Convert to Cartesian coordinates
    Ex = Eθ * cosϕ - Eϕ * sinϕ
    Ey = Eθ * sinϕ + Eϕ * cosϕ
    
    return Ex, Ey
end

"""
    calculate_fresnel_coefficients(kr2::Real, λ::Real, n_medium::Real, n_immersion::Real)

Calculates Fresnel transmission coefficients for light passing through an interface.

# Physics
- Handles dielectric interfaces between media with indices n₁ and n₂
- Uses Snell's law: n₁sinθ₁ = n₂sinθ₂
- Computes both p-(TM) and s-(TE) polarization coefficients
- Valid for arbitrary incidence angles below critical angle

# Parameters
- kr2: Radial spatial frequency squared
- λ: Wavelength in microns
- n_medium: Refractive index of medium
- n_immersion: Refractive index of immersion medium

# Returns
- FresnelP: p-polarization transmission coefficient
- FresnelS: s-polarization transmission coefficient
- sinθ₁: Sine of incident angle
- cosθ₁: Cosine of incident angle
"""
function calculate_fresnel_coefficients(kr2::Real, λ::Real, n_medium::Real, n_immersion::Real)
    # Calculate angles using complex sqrt for automatic handling
    sinθ₁ = sqrt(kr2)*λ/n_medium
    cosθ₁ = sqrt(complex(1 - kr2*λ^2/n_medium^2))
    cosθᵢ = sqrt(complex(1 - kr2*λ^2/n_immersion^2))

    # Calculate Fresnel coefficients
    FresnelP = 2.0*n_medium*cosθ₁/(n_medium*cosθᵢ + n_immersion*cosθ₁)
    FresnelS = 2.0*n_medium*cosθ₁/(n_medium*cosθ₁ + n_immersion*cosθᵢ)

    return FresnelP, FresnelS, sinθ₁, cosθ₁
end

"""
    Vector3DPupilPSF(nₐ::Real, λ::Real, pupil::PupilFunction; kwargs...)

Construct Vector3DPSF with pupil function calculated from vectorial diffraction theory.

# Arguments
- `nₐ`: Numerical aperture 
- `λ`: Wavelength in microns
- `pupil`: Base pupil function containing aberrations

# Keyword Arguments 
- `n_medium::Real=1.33`: Sample medium refractive index (water)
- `n_coverslip::Real=1.52`: Cover slip refractive index (glass)
- `n_immersion::Real=1.52`: Immersion medium refractive index (oil)
- `Σ::Union{Matrix,Nothing}=nothing`: 2×2 OTF rescaling matrix. Identity if nothing.
"""
function Vector3DPupilPSF(nₐ::Real, λ::Real, pupil::PupilFunction;
                         n_medium::Real=1.33,
                         n_coverslip::Real=1.52, 
                         n_immersion::Real=1.52,
                         Σ::Union{Matrix,Nothing}=nothing)
    
    grid_size = size(pupil.field, 1)
    xs = ys = range(-1, 1, length=grid_size)
    kmax = nₐ/λ
    kpix = 2kmax/(grid_size-1)
    
    # Initialize 6 pupil functions for field components
    pupils = Vector{PupilFunction}(undef, 6)
    for i in 1:6
        pupils[i] = PupilFunction(nₐ, λ, n_medium, 
                                zeros(Complex{Float64}, grid_size, grid_size))
    end
    
    # Fill pupil functions
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        
        if ρ > 1
            continue
        end
        
        # Get base aberrated field
        base_field = pupil.field[j,i]
        
        # Calculate angles
        ϕ = atan(y, x)
        
        # Before angle calculations, replace with:
        kr2 = (ρ * nₐ/λ)^2
        Tp, Ts, sinθ₁, cosθ₁ = calculate_fresnel_coefficients(kr2, λ, n_medium, n_immersion)
        
        # Calculate field components for each dipole orientation
        for (d, dipole) in enumerate([:x, :y, :z])
            # Calculate pupil field for this dipole
            Ex, Ey = calculate_pupil_field(ϕ, Tp, Ts, sinθ₁, cosθ₁;
                                         dvec=dipole==:x ? [1,0,0] :
                                              dipole==:y ? [0,1,0] : [0,0,1])
            
            # Combine with base aberrated field
            pupils[d].field[j,i] = Ex * base_field
            pupils[d+3].field[j,i] = Ey * base_field
        end
    end
    
    # Normalize pupils
    for p in pupils
        normalize!(p)
    end
    
    # Create VectorPupilFunction to handle the 6 components
    vpupil = VectorPupilFunction(nₐ, λ, n_medium, pupils)
    
    # Use default identity matrix for Σ if not provided
    Σ_mat = isnothing(Σ) ? Matrix{Float64}(I, 2, 2) : Σ
    
    # Create final Vector3DPupilPSF using main constructor
    Vector3DPupilPSF(nₐ, λ, n_medium, n_coverslip, n_immersion, Σ_mat, vpupil)
end

"""
    amplitude(p::Vector3DPupilPSF, x::Real, y::Real, z::Real, dipole::DipoleVector)

Calculate complex vector amplitude at given position for specific dipole orientation.
"""
function amplitude(p::Vector3DPupilPSF{T}, x::Real, y::Real, z::Real, dipole::DipoleVector) where T
    # Compute the pupil function field at the given position
    Ex, Ey = amplitude(p.pupil, x, y, z, dipole)
    return [Ex, Ey]
end

"""
    (psf::Vector3DPupilPSF)(x::Real, y::Real, z::Real, dipole::DipoleVector)

Evaluate PSF intensity for specific dipole orientation.
"""
function (psf::Vector3DPupilPSF)(x::Real, y::Real, z::Real, dipole::DipoleVector)
    E = amplitude(psf, x, y, z, dipole)
    return abs2(E[1]) + abs2(E[2])
end

"""
    (psf::Vector3DPupilPSF)(x::Real, y::Real, z::Real)

Evaluate PSF intensity for randomly oriented dipoles by summing over three orthogonal orientations.
"""
function (psf::Vector3DPupilPSF)(x::Real, y::Real, z::Real)
    # Calculate field for each orthogonal dipole orientation
    Ex = amplitude(psf, x, y, z, DipoleVector(1,0,0))
    Ey = amplitude(psf, x, y, z, DipoleVector(0,1,0))
    Ez = amplitude(psf, x, y, z, DipoleVector(0,0,1))
    
    # Sum intensities and normalize by number of orientations
    return (abs2(Ex[1]) + abs2(Ex[2]) + 
            abs2(Ey[1]) + abs2(Ey[2]) + 
            abs2(Ez[1]) + abs2(Ez[2])) / 3
end

"""
    integrate_pixels(psf::Vector3DPupilPSF, camera::AbstractCamera, 
                    emitter::DipoleEmitter3D; sampling::Integer=2)

Calculate integrated pixel values for vector PSF and dipole emitter.
"""
function integrate_pixels(psf::Vector3DPupilPSF, camera::AbstractCamera,
                        emitter::DipoleEmitter3D; sampling::Integer=2)
    return _integrate_pixels_generic(
        psf, camera, emitter,
        (p,x,y) -> p(x, y, emitter.z, emitter.dipole),
        typeof(emitter.photons);
        sampling=sampling
    )
end

"""
    integrate_pixels(psf::Vector3DPupilPSF, camera::AbstractCamera,
                    emitter::Emitter3D; sampling::Integer=2)

Calculate integrated pixel values for vector PSF and non-dipole emitter.
"""
function integrate_pixels(psf::Vector3DPupilPSF, camera::AbstractCamera,
                        emitter::Emitter3D; sampling::Integer=2)
    return _integrate_pixels_generic(
        psf, camera, emitter,
        (p,x,y) -> p(x, y, emitter.z),
        typeof(emitter.photons);
        sampling=sampling
    )
end

