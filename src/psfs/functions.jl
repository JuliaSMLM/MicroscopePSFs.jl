"""
    (psf::AbstractPSF)(x, y)

Evaluate normalized PSF intensity at position (x,y) relative to PSF center.

# Arguments
- `x, y`: Position in microns relative to PSF center

# Returns
- Intensity normalized to integrate to 1

# Coordinates
Input coordinates (x,y) are in physical units (microns) relative to PSF center.
"""
function (psf::AbstractPSF)(x::Real, y::Real)
    error("PSF evaluation not implemented for $(typeof(psf))")
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
"""
function amplitude(psf::AbstractPSF, x::Real, y::Real)
    error("Amplitude calculation not implemented for $(typeof(psf))")
end
