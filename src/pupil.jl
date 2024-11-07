"""
    PupilFunction{T}

Represents a pupil function with physical parameters.

# Fields
- `nₐ`: Numerical aperture
- `λ`: Wavelength in μm
- `n`: Refractive index
- `data`: Complex-valued pupil function array
"""
struct PupilFunction{T<:AbstractFloat}
    nₐ::T
    λ::T
    n::T
    data::Array{Complex{T},2}

    function PupilFunction(nₐ::Real, λ::Real, n::Real, data::Array{<:Complex})
        nₐ > 0 || throw(ArgumentError("Numerical aperture must be positive"))
        λ > 0 || throw(ArgumentError("Wavelength must be positive"))
        n > 0 || throw(ArgumentError("Refractive index must be positive"))
        T = promote_type(typeof(nₐ), typeof(λ), typeof(n), real(eltype(data)))
        new{T}(T(nₐ), T(λ), T(n), convert(Array{Complex{T}}, data))
    end
end

# Computed properties
"""Get maximum spatial frequency in μm⁻¹"""
kmax(p::PupilFunction) = p.nₐ / p.λ

"""Get central wavevector magnitude in μm⁻¹"""
k₀(p::PupilFunction) = p.n / p.λ

"""Get pupil plane sampling in μm⁻¹"""
kpixelsize(p::PupilFunction) = 2kmax(p) / (size(p.data,1) - 1)

# Constructor from Zernike coefficients
"""
    PupilFunction(nₐ::Real, λ::Real, n::Real, 
                 coeffs::AbstractVector, 
                 orders::AbstractVector{<:Tuple{Int,Int}};
                 grid_size::Int=64)

Create a PupilFunction from Zernike coefficients.
"""
function PupilFunction(nₐ::Real, λ::Real, n::Real,
                      coeffs::AbstractVector,
                      orders::AbstractVector{<:Tuple{Int,Int}};
                      grid_size::Int=64)
    T = promote_type(typeof(nₐ), typeof(λ), typeof(n), eltype(coeffs))
    data = zeros(Complex{T}, grid_size, grid_size)
    
    # Generate normalized coordinate grid
    xs = ys = range(-1, 1, length=grid_size)
    
    for i in 1:grid_size, j in 1:grid_size
        x, y = xs[i], ys[j]
        ρ = sqrt(x^2 + y^2)
        ρ > 1 && continue
        θ = atan(y, x)
        
        for (c, (n,m)) in zip(coeffs, orders)
            data[j,i] += c * zernike(n, m, ρ, θ)
        end
    end
    
    return PupilFunction(nₐ, λ, n, data)
end

# Normalization
"""
    normalize!(p::PupilFunction)

Normalize pupil function to unit energy using Parseval's theorem.
Returns the normalized PupilFunction.
"""
function normalize!(p::PupilFunction)
    # Energy normalization using Parseval's theorem
    total_energy = sum(abs2, p.data) * kpixelsize(p)^2
    p.data ./= sqrt(total_energy)
    return p
end

# Display and visualization
function Base.show(io::IO, p::PupilFunction)
    sz = size(p.data, 1)
    print(io, "PupilFunction(NA=$(p.nₐ), λ=$(p.λ)μm, n=$(p.n), $(sz)x$(sz))")
end

# Utility functions for pupil manipulation
"""
    apply_defocus!(p::PupilFunction, z::Real)

Apply defocus phase to pupil function for propagation distance z.
"""
function apply_defocus!(p::PupilFunction, z::Real)
    sz = size(p.data, 1)
    kpix = kpixelsize(p)
    k0_center = (sz + 1) ÷ 2
    
    for i in 1:sz, j in 1:sz
        kx = kpix * (i - k0_center)
        ky = kpix * (j - k0_center)
        kr2 = kx^2 + ky^2
        
        if kr2 < kmax(p)^2
            kz = sqrt(complex(k₀(p)^2 - kr2))
            p.data[j,i] *= exp(2π * im * z * kz)
        end
    end
    return p
end

"""
    apply_aperture!(p::PupilFunction, radius::Real=1.0)

Apply circular aperture to pupil function. Radius is relative to NA.
"""
function apply_aperture!(p::PupilFunction, radius::Real=1.0)
    sz = size(p.data, 1)
    xs = ys = range(-1, 1, length=sz)
    
    for i in 1:sz, j in 1:sz
        x, y = xs[i], ys[j]
        if sqrt(x^2 + y^2) > radius
            p.data[j,i] = 0
        end
    end
    return p
end

"""
Calculate complex amplitude from pupil function integration.

# Arguments
- `p::PupilFunction`: Pupil function
- `x::Real`: X position in μm
- `y::Real`: Y position in μm
- `z::Real`: Z position in μm

Returns complex amplitude at specified position.
"""
function amplitude(p::PupilFunction{T}, x::Real, y::Real, z::Real) where T
    sz = size(p.data, 1)
    kpix = kpixelsize(p)
    k0_center = (sz + 1) ÷ 2
    
    result = zero(Complex{T})
    kmax² = kmax(p)^2
    k₀² = k₀(p)^2
    
    @inbounds for i in 1:sz, j in 1:sz
        kx = kpix * (i - k0_center)
        ky = kpix * (j - k0_center)
        kr2 = kx^2 + ky^2
        
        if kr2 < kmax²
            kz = sqrt(complex(k₀² - kr2))
            phase = 2π * (x*kx + y*ky + z*kz)
            result += p.data[j,i] * exp(im * phase)
        end
    end
    
    return result * kpix^2
end

