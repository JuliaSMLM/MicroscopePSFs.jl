# src/vector_pupil.jl

"""
    VectorPupilFunction{T}

Vector pupil function storing Ex,Ey field components as PupilFunctions.

# Fields
- `nₐ::T`: Numerical aperture
- `λ::T`: Wavelength in μm
- `n_medium::T`: Refractive index of the sample medium
- `n_coverslip::T`: Refractive index of the coverslip
- `n_immersion::T`: Refractive index of the immersion medium
- `Ex::PupilFunction{T}`: x-component of electric field
- `Ey::PupilFunction{T}`: y-component of electric field
"""
struct VectorPupilFunction{T<:AbstractFloat}
    nₐ::T
    λ::T
    n_medium::T
    n_coverslip::T
    n_immersion::T
    Ex::PupilFunction{T}
    Ey::PupilFunction{T}

    function VectorPupilFunction(nₐ::Real, λ::Real, n_medium::Real, n_coverslip::Real, n_immersion::Real,
                               Ex::PupilFunction, Ey::PupilFunction)
        # Validation
        nₐ > 0 || throw(ArgumentError("Numerical aperture must be positive"))
        λ > 0 || throw(ArgumentError("Wavelength must be positive"))
        n_medium > 0 || throw(ArgumentError("n_medium must be positive"))
        n_coverslip > 0 || throw(ArgumentError("n_coverslip must be positive"))
        n_immersion > 0 || throw(ArgumentError("n_immersion must be positive"))
        
        size(Ex.field) == size(Ey.field) || 
            throw(ArgumentError("Ex and Ey pupils must have same size"))
        
        T = promote_type(typeof(nₐ), typeof(λ),
                        typeof(n_medium), typeof(n_coverslip), typeof(n_immersion),
                        real(eltype(Ex.field)), real(eltype(Ey.field)))
        
        new{T}(T(nₐ), T(λ), T(n_medium), T(n_coverslip), T(n_immersion), Ex, Ey)
    end
end

"""
    VectorPupilFunction(nₐ::Real, λ::Real, n_medium::Real, n_coverslip::Real, n_immersion::Real, grid_size::Integer)

Create a vector pupil function with empty field components.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in μm
- `n_medium`: Sample medium refractive index
- `n_coverslip`: Cover slip refractive index
- `n_immersion`: Immersion medium refractive index
- `grid_size`: Size of pupil grid (grid_size × grid_size)

# Returns
- `VectorPupilFunction` with zero-initialized field components
"""
function VectorPupilFunction(nₐ::Real, λ::Real, n_medium::Real, n_coverslip::Real, n_immersion::Real, grid_size::Integer)
    Ex = PupilFunction(nₐ, λ, n_medium, zeros(Complex{Float64}, grid_size, grid_size))
    Ey = PupilFunction(nₐ, λ, n_medium, zeros(Complex{Float64}, grid_size, grid_size))
    VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, Ex, Ey)
end

# Computed properties
"""Get maximum spatial frequency in μm⁻¹"""
kmax(p::VectorPupilFunction) = kmax(p.Ex)

"""Get central wavevector magnitude in μm⁻¹"""
k₀(p::VectorPupilFunction) = k₀(p.Ex)

"""Get pupil plane sampling in μm⁻¹"""
kpixelsize(p::VectorPupilFunction) = kpixelsize(p.Ex)

"""
    normalize!(p::VectorPupilFunction)

Normalize the electric field components of the vector pupil function.
Ensures total energy across both components equals 1.

# Arguments
- `p`: Vector pupil function to normalize

# Returns
- Normalized vector pupil function
"""
function normalize!(p::VectorPupilFunction)
    kpix² = kpixelsize(p)^2
    total_energy = (sum(abs2, p.Ex.field) + sum(abs2, p.Ey.field)) * kpix²
    scale = 1/sqrt(total_energy)
    
    p.Ex.field .*= scale
    p.Ey.field .*= scale
    return p
end

"""
    fill_vector_pupils!(vpupil::VectorPupilFunction, dipole::DipoleVector,
                      base_pupil::Union{Nothing, PupilFunction}=nothing;
                      normalize::Bool=true)

Fill vector pupil function with field components including dipole orientation,
base aberrations, and proper apodization.

This pre-calculates all position-independent factors of the pupil function.

# Arguments
- `vpupil`: Vector pupil function to fill
- `dipole`: Dipole orientation vector
- `base_pupil`: Optional base aberration pupil function

# Keyword Arguments
- `normalize`: Whether to normalize the pupil function (default: true)

# Returns
- Filled vector pupil function (normalized if `normalize=true`)
"""
function fill_vector_pupils!(vpupil::VectorPupilFunction, 
                           dipole::DipoleVector, 
                           base_pupil::Union{Nothing, PupilFunction}=nothing;
                           normalize::Bool=true)
    
    grid_size = size(vpupil.Ex.field, 1)
    xs = ys = range(-1, 1, length=grid_size)
    kmax = vpupil.nₐ/vpupil.λ
    
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        
        # Convert normalized pupil coordinates to k-space
        kr2 = (ρ * kmax)^2
        
        # Calculate angles and convert to ComplexF64 to handle evanescent waves
        ϕ = ComplexF64(atan(y, x))
            
        # Calculate angles in medium for field components
        sinθ = sqrt(complex(kr2 * vpupil.λ^2 / vpupil.n_medium^2))
        cosθ = sqrt(complex(1 - kr2 * vpupil.λ^2 / vpupil.n_medium^2))
        
        # Calculate Fresnel coefficients for interfaces
        Tp, Ts = calculate_fresnel_coefficients(
            kr2, vpupil.λ, vpupil.n_medium, vpupil.n_coverslip, vpupil.n_immersion)
        
        # Calculate field components for dipole orientation
        Ex, Ey = calculate_dipole_field_components(ϕ, sinθ, cosθ, dipole, Tp, Ts)
        
        # Calculate apodization factor
        apod = calculate_apodization(kr2, vpupil.λ, vpupil.n_medium, vpupil.n_immersion)
        
        # Apply apodization 
        Ex *= apod
        Ey *= apod
        
        # Combine with base aberration if provided
        if !isnothing(base_pupil)
            base_field = base_pupil.field[j,i]
            Ex *= base_field
            Ey *= base_field
        end
        
        # Store values in pupil
        vpupil.Ex.field[j,i] = Ex
        vpupil.Ey.field[j,i] = Ey
    end
    
    # Normalize pupils only if requested
    if normalize
        normalize!(vpupil)
    end
    
    return vpupil
end



# Display methods
function Base.show(io::IO, p::VectorPupilFunction)
    sz = size(p.Ex.field, 1)
    print(io, "VectorPupilFunction(NA=$(p.nₐ), λ=$(p.λ)μm, ",
          "n_medium=$(p.n_medium), n_coverslip=$(p.n_coverslip), ",
          "n_immersion=$(p.n_immersion), $(sz)x$(sz))")
end

"""
    copy(vp::VectorPupilFunction)

Return a deep copy of `vp` with independent `Ex` and `Ey` field arrays.
"""
Base.copy(vp::VectorPupilFunction) = VectorPupilFunction(
    vp.nₐ, vp.λ, vp.n_medium, vp.n_coverslip, vp.n_immersion,
    copy(vp.Ex), copy(vp.Ey))

"""
    apply_defocus!(vp::VectorPupilFunction, z::Real, z_stage::Real=zero(z))

Apply the vectorial axial defocus phase to a `VectorPupilFunction` in place.

The phase `exp(im·2π·(k_{z,medium}·z − k_{z,immersion}·z_stage))` is multiplied into
both the `Ex` and `Ey` field components at every in-aperture sample, using the optical
parameters (`nₐ`, `λ`, `n_medium`, `n_coverslip`, `n_immersion`) stored in `vp`. Samples
outside the aperture are left unchanged (they are already zero).

This reuses the same k-space grid and phase convention (`calculate_wave_vectors`,
`calculate_axial_phase`) as `VectorPSF` evaluation, so the result is the pupil whose
Fourier transform yields the field at emitter depth `z` and stage position `z_stage`.

# Arguments
- `vp`: Vector pupil function to modify
- `z`: Emitter depth above the coverslip (μm)
- `z_stage`: Distance the stage was moved from the nominal focal plane (μm, default 0)

# Returns
- The modified `vp`
"""
function apply_defocus!(vp::VectorPupilFunction, z::Real, z_stage::Real=zero(z))
    grid_size = size(vp.Ex.field, 1)
    kmax_val = vp.nₐ / vp.λ
    kpixel = 2 * kmax_val / (grid_size - 1)
    center = (grid_size + 1) / 2

    for i in 1:grid_size, j in 1:grid_size
        kx = (i - center) * kpixel
        ky = (j - center) * kpixel
        kr2 = kx^2 + ky^2
        kr2 > kmax_val^2 && continue

        kz_medium, kz_coverslip, kz_immersion = calculate_wave_vectors(
            kr2, vp.λ, vp.n_medium, vp.n_coverslip, vp.n_immersion)
        axial_phase = calculate_axial_phase(z, z_stage, kz_medium, kz_coverslip, kz_immersion)
        phase_factor = exp(im * 2π * axial_phase)

        vp.Ex.field[j, i] *= phase_factor
        vp.Ey.field[j, i] *= phase_factor
    end

    return vp
end
