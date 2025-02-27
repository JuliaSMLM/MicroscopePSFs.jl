# src/psfs/spline_psf.jl

"""
    SplinePSF{T<:AbstractFloat} <: AbstractPSF

PSF representation using cubic spline interpolation.
Provides fast, accurate evaluation of PSFs from measured or simulated data.

# Fields
- `coefficients::Array{T,7}`: Cubic spline coefficients [nx-1, ny-1, nz-1, 4, 4, 4]
- `x_knots::Vector{T}`: x-coordinate knot positions
- `y_knots::Vector{T}`: y-coordinate knot positions
- `z_knots::Vector{T}`: z-coordinate knot positions
"""
struct SplinePSF{T<:AbstractFloat} <: AbstractPSF
    coefficients::Array{T,7}
    x_knots::Vector{T}
    y_knots::Vector{T}
    z_knots::Vector{T}
end

"""
    SplinePSF(psf_stack::AbstractArray{<:Real,3};
              x_coords::AbstractVector,
              y_coords::AbstractVector,
              z_coords::AbstractVector)

Create a SplinePSF from a 3D PSF measurement stack.

# Arguments
- `psf_stack`: 3D array containing PSF intensity values [y, x, z]
- `x_coords`, `y_coords`, `z_coords`: Physical coordinates for the grid points

# Returns
- `SplinePSF` with precomputed cubic spline coefficients
"""
function SplinePSF(psf_stack::AbstractArray{<:Real,3};
                   x_coords::AbstractVector,
                   y_coords::AbstractVector,
                   z_coords::AbstractVector)
    
    # Normalize PSF to sum to 1
    normalized_psf = psf_stack / sum(psf_stack)
    
    # Calculate spline coefficients directly from the data
    T = eltype(normalized_psf)
    coefficients = calculate_spline_coefficients(normalized_psf, x_coords, y_coords, z_coords)
    
    return SplinePSF{T}(coefficients, x_coords, y_coords, z_coords)
end

"""
    SplinePSF(psf::AbstractPSF; 
              x_range::AbstractRange=range(-1.0, 1.0, length=21),
              y_range::AbstractRange=range(-1.0, 1.0, length=21),
              z_range::AbstractRange=range(-1.0, 1.0, length=21))

Create a SplinePSF from any other PSF type by direct sampling.

# Arguments
- `psf`: Source PSF to sample
- `x_range`, `y_range`, `z_range`: Coordinate ranges to sample

# Returns
- `SplinePSF` that evaluates much faster than the original PSF
"""
function SplinePSF(psf::AbstractPSF; 
                  x_range::AbstractRange=range(-1.0, 1.0, length=21),
                  y_range::AbstractRange=range(-1.0, 1.0, length=21),
                  z_range::AbstractRange=range(-1.0, 1.0, length=21))
    
    # Sample PSF directly at the desired knot positions
    psf_stack = Array{Float64}(undef, length(y_range), length(x_range), length(z_range))
    for (iz, z) in enumerate(z_range)
        for (iy, y) in enumerate(y_range)
            for (ix, x) in enumerate(x_range)
                psf_stack[iy, ix, iz] = psf(x, y, z)
            end
        end
    end
    
    # Create SplinePSF directly from sampled values
    return SplinePSF(psf_stack;
                  x_coords=collect(x_range),
                  y_coords=collect(y_range),
                  z_coords=collect(z_range))
end

# Implement AbstractPSF interface

"""
    (psf::SplinePSF)(x::Real, y::Real, z::Real)

Evaluate the spline PSF at a 3D position using the precomputed coefficients.
"""
function (psf::SplinePSF)(x::Real, y::Real, z::Real)
    # Find the grid cell containing (x,y,z)
    ix, iy, iz = find_cell(psf, x, y, z)
    
    # Return zero if outside the domain
    if isnothing(ix)
        return zero(eltype(psf.coefficients))
    end
    
    # Get local coordinates within the cell (normalized to [0,1])
    x_norm = (x - psf.x_knots[ix]) / (psf.x_knots[ix+1] - psf.x_knots[ix])
    y_norm = (y - psf.y_knots[iy]) / (psf.y_knots[iy+1] - psf.y_knots[iy])
    z_norm = (z - psf.z_knots[iz]) / (psf.z_knots[iz+1] - psf.z_knots[iz])
    
    # Evaluate the cubic polynomial using the stored coefficients
    return evaluate_tricubic(psf.coefficients[ix, iy, iz, :, :, :], x_norm, y_norm, z_norm)
end

"""
    (psf::SplinePSF)(x::Real, y::Real)

Evaluate the spline PSF at a 2D position in the focal plane (z=0).
"""
function (psf::SplinePSF)(x::Real, y::Real)
    return psf(x, y, 0.0)
end

"""
    amplitude(psf::SplinePSF, x::Real, y::Real, z::Real)

Calculate amplitude at position (x,y,z).
For real-valued PSFs, this is the square root of intensity.
"""
function amplitude(psf::SplinePSF, x::Real, y::Real, z::Real)
    intensity = psf(x, y, z)
    return sqrt(complex(intensity))  # Use complex to handle negative values that might occur from interpolation errors
end

"""
    amplitude(psf::SplinePSF, x::Real, y::Real)

Calculate amplitude at position (x,y) in focal plane.
"""
function amplitude(psf::SplinePSF, x::Real, y::Real)
    return amplitude(psf, x, y, 0.0)
end

"""
    find_cell(psf::SplinePSF, x::Real, y::Real, z::Real)

Find the grid cell indices containing position (x,y,z).
Returns (ix, iy, iz) or (nothing, nothing, nothing) if outside domain.
"""
function find_cell(psf::SplinePSF, x::Real, y::Real, z::Real)
    # Check if position is within domain
    if x < minimum(psf.x_knots) || x >= maximum(psf.x_knots) ||
       y < minimum(psf.y_knots) || y >= maximum(psf.y_knots) ||
       z < minimum(psf.z_knots) || z >= maximum(psf.z_knots)
        return nothing, nothing, nothing
    end
    
    # Find indices using binary search
    ix = max(1, searchsortedlast(psf.x_knots, x))
    iy = max(1, searchsortedlast(psf.y_knots, y))
    iz = max(1, searchsortedlast(psf.z_knots, z))
    
    # Handle edge case - if at last knot
    if ix == length(psf.x_knots) ix -= 1 end
    if iy == length(psf.y_knots) iy -= 1 end
    if iz == length(psf.z_knots) iz -= 1 end
    
    return ix, iy, iz
end

# IO methods

"""
    save_spline_psf(filename::String, psf::SplinePSF)

Save a SplinePSF to HDF5 file.
"""
function save_spline_psf(filename::String, psf::SplinePSF)
    h5open(filename, "w") do file
        file["coefficients"] = psf.coefficients
        file["x_knots"] = psf.x_knots
        file["y_knots"] = psf.y_knots
        file["z_knots"] = psf.z_knots
    end
end

"""
    load_spline_psf(filename::String)

Load a SplinePSF from HDF5 file.
"""
function load_spline_psf(filename::String)
    h5open(filename, "r") do file
        coefficients = read(file["coefficients"])
        x_knots = read(file["x_knots"])
        y_knots = read(file["y_knots"])
        z_knots = read(file["z_knots"])
        
        T = eltype(coefficients)
        return SplinePSF{T}(coefficients, x_knots, y_knots, z_knots)
    end
end

# Pretty printing
function Base.show(io::IO, psf::SplinePSF)
    nx = length(psf.x_knots) - 1
    ny = length(psf.y_knots) - 1
    nz = length(psf.z_knots) - 1
    
    x_range = maximum(psf.x_knots) - minimum(psf.x_knots)
    z_range = maximum(psf.z_knots) - minimum(psf.z_knots)
    
    print(io, "SplinePSF($(nx)×$(ny)×$(nz) grid, $(round(x_range, digits=2))×$(round(z_range, digits=2))μm)")
end