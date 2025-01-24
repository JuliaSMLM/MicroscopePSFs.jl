# src/pupil.jl

"""
    PupilFunction{T}

Represents a pupil function with physical parameters.

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in μm
- `n`: Refractive index
- `field`: Complex-valued pupil function array
"""
struct PupilFunction{T<:AbstractFloat}
    nₐ::T
    λ::T
    n::T
    field::Array{Complex{T},2}

    function PupilFunction(nₐ::Real, λ::Real, n::Real, field::Array{<:Complex})
        nₐ > 0 || throw(ArgumentError("Numerical aperture must be positive"))
        λ > 0 || throw(ArgumentError("Wavelength must be positive"))
        n > 0 || throw(ArgumentError("Refractive index must be positive"))
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), real(eltype(field)))
        new{T}(T(nₐ), T(λ), T(n), convert(Array{Complex{T}}, field))
    end
end

# Computed properties
"""
    kmax(p::PupilFunction)

Get maximum spatial frequency in μm⁻¹"""
kmax(p::PupilFunction) = p.nₐ / p.λ

"""
    k₀(p::PupilFunction)

Get central wavevector magnitude in μm⁻¹"""
k₀(p::PupilFunction) = p.n / p.λ

"""
    kpixelsize(p::PupilFunction)

Get pupil plane sampling in μm⁻¹"""
kpixelsize(p::PupilFunction) = 2kmax(p) / (size(p.field, 1) - 1)

"""
    PupilFunction(nₐ::Real, λ::Real, n::Real, 
                  zc::ZernikeCoefficients;
                  grid_size::Int=64)

Create a PupilFunction from ZernikeCoefficients type.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `n`: Refractive index
- `zc`: ZernikeCoefficients containing amplitude and phase coefficients

# Keyword Arguments
- `grid_size`: Size of square sampling grid (default: 64)

# Returns
- `PupilFunction` with complex field computed from Zernike polynomials
"""
function PupilFunction(nₐ::Real, λ::Real, n::Real,
    zc::ZernikeCoefficients;
    grid_size::Int=64)
    T = promote_type(typeof(nₐ), typeof(λ), typeof(n), eltype(zc.phase))
    field = zeros(Complex{T}, grid_size, grid_size)

    # Generate normalized coordinate grid
    xs = ys = range(-1, 1, length=grid_size)

    # Maximum radial order from coefficient length
    max_n = max_radial_order(length(zc.phase))

    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        if ρ > 1
            continue
        end
        θ = atan(y, x)

        # Initialize amplitude and phase at this point
        amplitude = 1.0
        total_phase = 0.0

        # Accumulate amplitude aberrations
        for n_rad in 0:max_n
            for m in -n_rad:2:n_rad
                idx = nl2osa(n_rad, m) + 1  # +1 for 1-based indexing
                if idx > length(zc.mag)
                    continue
                end

                # Retrieve the amplitude coefficient
                amp_coeff = zc.mag[idx]

                # Skip negligible coefficients
                if abs(amp_coeff) < 1e-10
                    continue
                end

                # Evaluate Zernike polynomial
                Z = zernikepolynomial(n_rad, m, ρ, θ)

                # Accumulate amplitude contribution
                amplitude += amp_coeff * Z
            end
        end

        # Accumulate phase aberrations
        for n_rad in 0:max_n
            for m in -n_rad:2:n_rad
                idx = nl2osa(n_rad, m) + 1
                if idx > length(zc.phase)
                    continue
                end

                # Retrieve the phase coefficient (in radians)
                phase_coeff = zc.phase[idx]

                # Skip negligible coefficients
                if abs(phase_coeff) < 1e-10
                    continue
                end

                # Evaluate Zernike polynomial
                Z = zernikepolynomial(n_rad, m, ρ, θ)

                # Accumulate phase contribution
                total_phase += phase_coeff * Z
            end
        end

        # Compute the complex field at this point
        field[j, i] = amplitude * exp(im * total_phase)
    end

    return PupilFunction(nₐ, λ, n, field)
end


# Normalization
"""
    normalize!(p::PupilFunction)

Normalize pupil function to unit energy using Parseval's theorem.
"""
function normalize!(p::PupilFunction)
    # Energy normalization using Parseval's theorem
    total_energy = sum(abs2, p.field) * kpixelsize(p)^2
    p.field ./= sqrt(total_energy)
    return p
end

# Display and visualization
function Base.show(io::IO, p::PupilFunction)
    sz = size(p.field, 1)
    print(io, "PupilFunction(NA=$(p.nₐ), λ=$(p.λ)μm, n=$(p.n), $(sz)x$(sz))")
end

# Utility functions for pupil manipulation
"""
    apply_defocus!(p::PupilFunction, z::Real)

Apply defocus phase to pupil function for propagation distance z.
"""
function apply_defocus!(p::PupilFunction, z::Real)
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    k0_center = (sz + 1) / 2

    for i in 1:sz, j in 1:sz
        kx = kpix * (i - k0_center)
        ky = kpix * (j - k0_center)
        kr2 = kx^2 + ky^2

        if kr2 < kmax(p)^2
            kz = sqrt(complex(k₀(p)^2 - kr2))
            p.field[j, i] *= exp(2π * im * z * kz)
        end
    end
    return p
end

"""
    apply_aperture!(p::PupilFunction, radius::Real=1.0)

Apply circular aperture to pupil function. Radius is relative to NA.
"""
function apply_aperture!(p::PupilFunction, radius::Real=1.0)
    sz = size(p.field, 1)
    xs = ys = range(-1, 1, length=sz)

    for i in 1:sz, j in 1:sz
        x, y = xs[i], ys[j]
        if sqrt(x^2 + y^2) > radius
            p.field[j, i] = 0
        end
    end
    return p
end

"""
    amplitude(p::PupilFunction, x::Real, y::Real, z::Real)

Calculate complex amplitude from pupil function integration.

# Arguments
- `p::PupilFunction`: Pupil function
- `x::Real`: X position in μm
- `y::Real`: Y position in μm
- `z::Real`: Z position in μm"""
function amplitude(p::PupilFunction{T}, x::Real, y::Real, z::Real) where {T}
    sz = size(p.field, 1)
    kpix = kpixelsize(p)
    # k0_center = (sz + 1) ÷ 2
    k0_center = (sz + 1) / 2

    result = zero(Complex{T})
    kmax² = kmax(p)^2
    k₀² = k₀(p)^2

    @inbounds for i in 1:sz, j in 1:sz
        kx = kpix * (i - k0_center)
        ky = kpix * (j - k0_center)
        kr2 = kx^2 + ky^2

        if kr2 < kmax²
            kz = sqrt(complex(k₀² - kr2))
            phase = 2π * (x * kx + y * ky + z * kz)
            result += p.field[j, i] * exp(im * phase)
        end
    end

    return result * kpix^2
end

"""
    VectorPupilFunction{T}

Vector pupil function with Ex,Ey components from x,y,z dipoles.

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in μm
- `n`: Refractive index
- `pupils`: Array of 6 PupilFunctions for [Ex_x,Ex_y,Ex_z,Ey_x,Ey_y,Ey_z]
"""
struct VectorPupilFunction{T<:AbstractFloat} 
    nₐ::T
    λ::T
    n::T
    pupils::Vector{PupilFunction{T}}

    function VectorPupilFunction(nₐ::Real, λ::Real, n::Real, 
                               pupils::Vector{<:PupilFunction})
        length(pupils) == 6 || throw(ArgumentError("Must provide 6 pupil functions"))
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), 
                        real(eltype(first(pupils).field)))
        new{T}(T(nₐ), T(λ), T(n), convert(Vector{PupilFunction{T}}, pupils))
    end
end


"""
    amplitude(p::VectorPupilFunction, x::Real, y::Real, z::Real; orientation=:x)
    
Calculate vector amplitude at given position
    
# Arguments
- `p::VectorPupilFunction`: Vector pupil function
- `x::Real`: X position in μm
- `y::Real`: Y position in μm 
- `z::Real`: Z position in μm

# Keyword Arguments
- `orientation`: Dipole orientation (:x, :y, or :z), defaults to :x
    
# Returns 
- `Vector{Complex}`: [Ex,Ey] components of E-field"""
function amplitude(p::VectorPupilFunction{T}, x::Real, y::Real, z::Real;
                  orientation::Symbol=:x) where {T}
    if orientation == :x
        Ex = amplitude(p.pupils[1], x, y, z)
        Ey = amplitude(p.pupils[4], x, y, z)
    elseif orientation == :y
        Ex = amplitude(p.pupils[2], x, y, z) 
        Ey = amplitude(p.pupils[5], x, y, z)
    elseif orientation == :z 
        Ex = amplitude(p.pupils[3], x, y, z)
        Ey = amplitude(p.pupils[6], x, y, z)
    else
        throw(ArgumentError("orientation must be :x, :y or :z"))
    end
    return [Ex, Ey]
end

