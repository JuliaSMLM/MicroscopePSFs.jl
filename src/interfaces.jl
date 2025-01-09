# src/interfaces.jl

"""
    (psf::AbstractPSF)(x, y)

Evaluate normalized PSF intensity at position (x,y) relative to PSF center.

# Arguments
- `x, y`: Position in microns relative to PSF center

# Returns
- Intensity normalized to integrate to 1

# Coordinates
Input coordinates (x,y) are in physical units (microns) relative to PSF center.

# Examples
```julia
psf = Airy2D(1.4, 0.532)  # NA=1.4, λ=532nm
intensity = psf(0.1, 0.2)  # Evaluate at (0.1μm, 0.2μm)
```
"""
function (psf::AbstractPSF)(x::Real, y::Real)
    error("PSF evaluation not implemented for $(typeof(psf))")
end

"""
    (psf::AbstractPSF)(x, y, z)

Evaluate PSF intensity at 3D position. z is axial distance from focus in microns.
"""
function (psf::AbstractPSF)(x::Real, y::Real, z::Real)
    error("3D PSF evaluation not implemented for $(typeof(psf))")
end

"""
    amplitude(psf::AbstractPSF, x::Real, y::Real)

Evaluate complex field amplitude at position (x,y) relative to PSF center.

# Arguments
- `x, y`: Position in microns relative to PSF center

# Returns
- Complex amplitude normalized such that |amplitude|² gives normalized intensity

# Coordinates
Input coordinates (x,y) are in physical units (microns) relative to PSF center.

# Examples
```julia
psf = Airy2D(1.4, 0.532)
amp = amplitude(psf, 0.1, 0.2)
intensity = abs2(amp)  # Convert to intensity
```
"""
function amplitude(psf::AbstractPSF, x::Real, y::Real)
    error("Amplitude calculation not implemented for $(typeof(psf))")
end

"""
    amplitude(psf::AbstractPSF, x::Real, y::Real, z::Real)

Evaluate complex field amplitude in 3D. z is axial distance from focus in microns.
"""
function amplitude(psf::AbstractPSF, x::Real, y::Real, z::Real)
    error("3D amplitude calculation not implemented for $(typeof(psf))")
end

"""
    integrate_pixels(
        psf::AbstractPSF, 
        camera::AbstractCamera, 
        emitter::AbstractEmitter;
        sampling::Integer=2
    )

Calculate integrated pixel values for given PSF, camera geometry, and emitter position. 
This is the main interface for simulating microscope images of point emitters.

# Arguments
- `psf`: Point spread function (e.g., Airy2D, Scalar3DPSF)
- `camera`: Camera geometry defining pixel edges in microns
- `emitter`: Emitter with position and intensity
- `sampling`: Subpixel sampling factor for integration accuracy (default: 2)

# Returns
- Array of pixel values with dimensions [y,x] (matrix indexed [row,col])
- Values are normalized to sum to 1
- Array indices start at [1,1] for top-left pixel

# Camera Coordinates
The camera uses a coordinate system where:
- (0,0) is at the top-left corner
- x increases to the right
- y increases downward
- All units are in microns

# Examples
```julia
# Setup microscope components
camera = IdealCamera(0:0.1:2.0, 0:0.1:2.0)  # 20x20 pixels, 100nm size
psf = Airy2D(1.4, 0.532)  # NA=1.4, λ=532nm
emitter = Emitter2D(1.0, 1.0, 1000.0)  # At (1μm, 1μm) with 1000 photons

# Generate image
pixels = integrate_pixels(psf, camera, emitter)
```

See also: [`AbstractPSF`](@ref), [`AbstractCamera`](@ref), [`AbstractEmitter`](@ref)
"""
function integrate_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2)
    error("Pixel integration not implemented for $(typeof(psf))")
end
