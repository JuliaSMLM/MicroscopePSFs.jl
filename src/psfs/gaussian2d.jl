# src/psfs/gaussian2d.jl

# Type defined in types.jl 

"""
    (psf::GaussianPSF)(x, y)

Gaussian PSF intensity. Radially symmetric with standard deviation σ.
"""
function (psf::GaussianPSF{T})(x::Real, y::Real) where {T}
    C = T === Float32 && x isa Float32 && y isa Float32 ? Float32 : Float64
    if C === Float32
        width² = Float64(psf.σ)^2
        radius² = Float64(x)^2 + Float64(y)^2
        return Float32(exp(-radius² / (2width²)) / (2π * width²))
    end
    r² = x^2 + y^2
    return 1 / (2C(π) * psf.σ^2) * exp(-r² / (2psf.σ^2))
end

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
    x_integrals = [
        erf(
            (pixel_edges_x[i] - emitter.x) / scale,
            (pixel_edges_x[i + 1] - emitter.x) / scale
        ) / 2 for i in 1:nx
    ]
    y_integrals = [
        erf(
            (pixel_edges_y[j] - emitter.y) / scale,
            (pixel_edges_y[j + 1] - emitter.y) / scale
        ) / 2 for j in 1:ny
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


function amplitude(psf::GaussianPSF, x::Real, y::Real)
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