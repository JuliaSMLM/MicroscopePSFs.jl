# Interface

MicroscopePSFs.jl provides a consistent interface across all PSF types, making it easy to switch between different models with minimal code changes.

## Core Interface

All PSF types implement the following core interface methods:

### PSF Evaluation

```julia
# Evaluate PSF at specified position (call operator)
(psf::AbstractPSF)(x, y)      # 2D evaluation
(psf::AbstractPSF)(x, y, z)   # 3D evaluation

# Get complex amplitude
amplitude(psf::AbstractPSF, x, y)      # 2D amplitude
amplitude(psf::AbstractPSF, x, y, z)   # 3D amplitude
```

The coordinate system uses physical units (microns) with the origin at the PSF center.

### Pixel Integration

```julia
# Create camera and emitter
camera = IdealCamera(20, 20, 0.1)  # 20x20 pixels, 100nm pixel size
emitter = Emitter2D(1.0, 1.0, 1000.0)  # At (1μm, 1μm) with 1000 photons

# Generate realistic microscope image
pixels = integrate_pixels(psf, camera, emitter)

# For complex amplitude integration
pixels_amplitude = integrate_pixels_amplitude(psf, camera, emitter)
```

The `integrate_pixels` function is the primary way to generate physically realistic microscope images, accounting for:
- PSF shape
- Camera pixel geometry
- Emitter position and intensity

### Optimization with Support Regions

For performance optimization, both `integrate_pixels` and `integrate_pixels_amplitude` accept an optional `support` parameter to limit calculations to regions where the PSF has significant contribution:

```julia
# Calculate only within a radius of 0.5μm around the emitter
pixels = integrate_pixels(psf, camera, emitter, support=0.5)

# Or specify an explicit region
region = (-1.0, 1.0, -1.0, 1.0)  # (x_min, x_max, y_min, y_max) in μm
pixels = integrate_pixels(psf, camera, emitter, support=region)
```

This optimization is particularly valuable for multi-emitter simulations.

### Multi-Emitter Simulation

All PSF types support simulating images with multiple emitters:

```julia
# Create multiple emitters
emitters = [
    Emitter2D(1.0, 1.0, 1000.0),
    Emitter2D(1.5, 1.2, 800.0)
]

# Generate camera image with all emitters
pixels = integrate_pixels(psf, camera, emitters)

# Use support region for efficiency
pixels = integrate_pixels(psf, camera, emitters, support=0.5)
```

For incoherent emitters (typical fluorescence), intensities are summed. For coherent simulations, use `integrate_pixels_amplitude`.

## Working with PSFs

### Creating PSF Instances

Each PSF type has its own constructor with parameters specific to that model. For example:

```julia
# Create a GaussianPSF with sigma=150nm
psf_gaussian = GaussianPSF(0.15)

# Create an AiryPSF with NA=1.4 and wavelength=532nm
psf_airy = AiryPSF(1.4, 0.532)

# Create a VectorPSF with more parameters
# Note: Create a dipole vector for the orientation (z-axis in this case)
dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Dipole along z-axis
psf_vector = VectorPSF(
    1.4,                # Numerical aperture
    0.68,               # Wavelength in microns
    dipole_z,           # Dipole orientation (along optical axis)
    n_medium=1.33,      # Sample medium refractive index
    n_immersion=1.518,  # Immersion medium refractive index
    n_coverslip=1.518   # Coverslip refractive index
)
```

### Using the Same Code with Different PSF Models

The common interface allows you to write generic code that works with any PSF type:

```julia
function analyze_psf_width(psf::AbstractPSF)
    # Generate intensity profile
    x = range(-1, 1, length=100)
    intensities = [psf(xi, 0.0) for xi in x]
    
    # Calculate FWHM or other properties
    # ...
    
    return results
end

# Works with any PSF model that implements the 2D interface
# Note: VectorPSF and ScalarPSF only support 3D evaluation with (x,y,z)
results_gaussian = analyze_psf_width(GaussianPSF(0.15))
results_airy = analyze_psf_width(AiryPSF(1.4, 0.532))
```

## Performance Considerations

The interface is designed to be both flexible and performant. For high-performance applications, consider:

1. Using the `support` parameter to limit calculations to relevant regions
2. Using SplinePSF to pre-compute complex PSFs for faster evaluation
3. Selecting the appropriate PSF model for your needs (simpler models are faster)
