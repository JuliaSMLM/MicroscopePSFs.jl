# src/psfs/scalar3d.jl

"""
    ScalarPSF(nₐ::Real, λ::Real, n::Real;
                pupil=nothing,
                pupil_data::Union{Nothing,AbstractMatrix}=nothing,
                zernike_coeffs=nothing,
                grid_size::Integer=128)

Create a scalar 3D PSF with the specified optical parameters.

# Arguments
- `nₐ`: Numerical aperture
- `λ`: Wavelength in microns
- `n`: Refractive index of the medium

# Keyword Arguments
- `pupil`: Optional pre-created PupilFunction
- `pupil_data`: Optional complex matrix to initialize the pupil function
- `zernike_coeffs`: Optional ZernikeCoefficients for aberrations
- `grid_size`: Generated pupil grid size (default: 128); a supplied `pupil` or
  `pupil_data` retains its own grid. The intensity repeats in x and y with period
  `L = (grid_size - 1) * λ / (2nₐ)` in microns, using the actual pupil grid size.
  Keep render `support` below `L/2` or the frame narrower than `L`. A larger grid
  widens `L`. For NA 1.4 and λ 0.6 microns, grid 128 gives `L ≈ 27.21` microns.
  Grid 128 takes about four times the evaluation time of grid 64. See the
  ScalarPSF guide for the optics, accuracy, and tabulation measurements.

# Returns
- `ScalarPSF`: 3D PSF with scalar diffraction model

# Notes
- Exactly one of `pupil`, `pupil_data`, or `zernike_coeffs` should be provided
- If none are provided, creates an unaberrated pupil
- The constructed PSF stores the provided Zernike coefficients for later use

# Examples
```julia
# Create an unaberrated PSF
psf = ScalarPSF(1.4, 0.532, 1.518)

# Create PSF with Zernike aberrations
zc = ZernikeCoefficients(15)
add_astigmatism!(zc, 0.5)  # Add astigmatism
psf = ScalarPSF(1.4, 0.532, 1.518; zernike_coeffs=zc)
```
"""
function ScalarPSF(
        nₐ::Real, λ::Real, n::Real;
        pupil::Union{Nothing, PupilFunction} = nothing,
        pupil_data::Union{Nothing, AbstractMatrix} = nothing,
        zernike_coeffs::Union{Nothing, ZernikeCoefficients} = nothing,
        grid_size::Integer = 128
    )
    
    # Add this check to enforce n ≥ NA
    n >= nₐ || throw(ArgumentError("Refractive index (n=$n) must be ≥ numerical aperture (NA=$nₐ). " * 
                                  "This is a physical constraint since NA = n·sin(θ). " *
                                  "Consider using VectorPSF for proper handling of index mismatch."))
    
    # Rest of the constructor remains the same
    pupil_func = if !isnothing(pupil)
        pupil
    elseif !isnothing(pupil_data)
        PupilFunction(nₐ, λ, n, pupil_data)
    else
        zernike_coeffs_out = isnothing(zernike_coeffs) ? ZernikeCoefficients(1) : zernike_coeffs
        PupilFunction(nₐ, λ, n, zernike_coeffs_out; grid_size = Int(grid_size))
    end

    normalize!(pupil_func)
    stored_coeffs = !isnothing(zernike_coeffs) ? zernike_coeffs : nothing
    
    return ScalarPSF(nₐ, λ, n, pupil_func, stored_coeffs)
end

"""
    amplitude(psf::ScalarPSF{T}, x::Real, y::Real, z::Real) where {T}

Calculate complex amplitude at a 3D position using scalar diffraction theory.

# Arguments
- `psf`: ScalarPSF instance
- `x, y, z`: Position in microns relative to PSF center

# Returns
- Complex amplitude at the specified position

# Notes
- Uses Fourier optics to propagate from pupil to image space
- Accounts for defocus and aberrations encoded in the pupil function
"""
function amplitude(psf::ScalarPSF{T}, x::Real, y::Real, z::Real) where {T}
    # Use existing PupilFunction amplitude method
    return amplitude(psf.pupil, x, y, z)
end

"""
    (psf::ScalarPSF)(x::Real, y::Real, z::Real)

Evaluate PSF intensity at the given 3D position.

# Arguments
- `x, y, z`: Position in microns relative to PSF center

# Returns
- Intensity value at the specified position

# Notes
- Calculated as |amplitude|² of the complex field
- Follows standard operator syntax: `intensity = psf(x, y, z)`
"""
function (psf::ScalarPSF)(x::Real, y::Real, z::Real)
    return abs2(amplitude(psf, x, y, z))
end

# Generic pixel integration handles z-coordinates for scalar PSFs.

dipole_fields(psf::ScalarPSF, x::Real, y::Real, z::Real) = [[amplitude(psf, x, y, z)]]

function Base.show(io::IO, psf::ScalarPSF)
    print(io, "ScalarPSF(NA=$(psf.nₐ), λ=$(psf.λ)μm, n=$(psf.n))")
    if !isnothing(psf.zernike_coeffs)
        print(io, " with $(length(psf.zernike_coeffs)) Zernike terms")
    end
end
