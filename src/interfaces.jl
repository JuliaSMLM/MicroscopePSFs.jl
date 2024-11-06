
"""
    (psf::AbstractPSF)(x, y)

Evaluate the normalized PSF at physical coordinates (x,y).
Returns intensity value normalized to integrate to 1.
"""
function (psf::AbstractPSF)(x::Real, y::Real)
    error("PSF evaluation not implemented for $(typeof(psf))")
end

"""
    amplitude(psf::AbstractPSF, x::Real, y::Real)

Evaluate the complex field amplitude at physical coordinates (x,y).
"""
function amplitude(psf::AbstractPSF, x::Real, y::Real)
    error("Amplitude calculation not implemented for $(typeof(psf))")
end

"""
    integrated_pixels(
        psf::AbstractPSF, 
        camera::AbstractCamera, 
        emitter::AbstractEmitter;
        sampling::Integer=2
    )

Calculate integrated pixel values for the PSF given camera geometry and emitter.
`sampling` controls subpixel sampling for integration accuracy.

Returns: Array of pixel values normalized to sum to 1 
"""
function integrated_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    sampling::Integer=2
)
    error("Pixel integration not implemented for $(typeof(psf))")
end
