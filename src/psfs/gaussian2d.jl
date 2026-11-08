# src/psfs/gaussian2d.jl

# Type defined in types.jl 

"""
    (psf::GaussianPSF)(x, y)

Gaussian PSF intensity. Radially symmetric with standard deviation σ.
"""
function (psf::GaussianPSF{T})(x::Real, y::Real) where {T}
    C = T === Float32 && x isa Float32 && y isa Float32 ? Float32 : Float64
    if C === Float32
        return Float32(_gaussian_intensity_float64(psf, x, y))
    end
    σ = _widen_for_square(psf.σ)
    r² = _widen_for_square(x)^2 + _widen_for_square(y)^2
    return 1 / (2C(π) * σ^2) * exp(-r² / (2σ^2))
end

function _gaussian_intensity_float64(psf::GaussianPSF{Float32}, x::Real, y::Real)
    width² = Float64(psf.σ)^2
    radius² = Float64(x)^2 + Float64(y)^2
    return exp(-radius² / (2width²)) / (2π * width²)
end

_amplitude_float64(psf::GaussianPSF{Float32}, x::Real, y::Real) =
    sqrt(_gaussian_intensity_float64(psf, x, y))

# All intensity integration routes dispatch here after choosing their pixel region.
function integrate_pixels!(
        result::AbstractMatrix{T},
        psf::GaussianPSF,
        pixel_edges_x::AbstractVector,
        pixel_edges_y::AbstractVector,
        emitter::AbstractEmitter;
        sampling::Integer = 2,
        threaded::Bool = true
    ) where {T <: Real}
    nx, ny = length(pixel_edges_x) - 1, length(pixel_edges_y) - 1
    size(result) == (ny, nx) || throw(
        DimensionMismatch(
            "Result array dimensions $(size(result)) do not match expected ($ny, $nx)"
        )
    )
    scale = sqrt(typeof(psf.σ)(2)) * psf.σ
    # Integer edges or positions become floating point, so an offset cannot wrap.
    edges_x, edges_y = _float_coordinates(pixel_edges_x), _float_coordinates(pixel_edges_y)
    x₀, y₀ = _float_coordinate(emitter.x), _float_coordinate(emitter.y)
    x_integrals = [
        erf((edges_x[i] - x₀) / scale, (edges_x[i + 1] - x₀) / scale) / 2 for i in 1:nx
    ]
    y_integrals = [
        erf((edges_y[j] - y₀) / scale, (edges_y[j + 1] - y₀) / scale) / 2 for j in 1:ny
    ]
    function integrate_column(i)
        for j in 1:ny
            result[j, i] = emitter.photons * x_integrals[i] * y_integrals[j]
        end
        return
    end
    if threaded
        Threads.@threads for i in 1:nx
            integrate_column(i)
        end
    else
        for i in 1:nx
            integrate_column(i)
        end
    end
    return result
end


function amplitude(psf::GaussianPSF{T}, x::Real, y::Real) where {T}
    if T === Float32 && x isa Float32 && y isa Float32
        return Float32(_amplitude_float64(psf, x, y))
    end
    return sqrt(psf(x, y))
end




"""
    GaussianPSF(psf::AiryPSF)

Convert an AiryPSF to an approximate GaussianPSF by matching the FWHM.

Uses the relationship that σ ≈ 0.22 * λ/NA
"""
function GaussianPSF(psf::AiryPSF)
    σ = 0.22 * psf.λ / psf.nₐ
    return GaussianPSF(σ)
end

Base.show(io::IO, psf::GaussianPSF) = print(io, "GaussianPSF(σ=$(psf.σ)μm)")