# Gaussian2D

The `Gaussian2D` PSF model represents the microscope point spread function as a simple isotropic 2D Gaussian function. While this is a mathematical approximation rather than a physical model, it works well for many applications and is computationally very efficient.

## Mathematical Model

The Gaussian2D PSF is defined as:

```math
I(x, y) = \frac{1}{2\pi\sigma^2} \exp\left(-\frac{x^2 + y^2}{2\sigma^2}\right)
```

where:
- ``x, y`` are coordinates in physical units (typically microns)
- ``\sigma`` is the standard deviation in the same units

This function is normalized to integrate to 1 over the entire domain.

## Constructor

```julia
Gaussian2D(σ)
```

### Parameters
- `σ`: Standard deviation in microns

### Examples

```julia
# Create a Gaussian PSF with 150nm standard deviation
psf = Gaussian2D(0.15)

# Create a Gaussian approximation of an Airy disk
airy_psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm
gaussian_approximation = Gaussian2D(airy_psf)  # Automatically sets appropriate σ
```

## Methods

### Evaluation

```julia
# Evaluate PSF at a specific position
intensity = psf(x, y)

# Get complex amplitude (returns sqrt of intensity for Gaussian)
amp = amplitude(psf, x, y)
```

### Creating Images

```julia
# Create a grid of positions
x = range(-1, 1, length=100)
y = range(-1, 1, length=100)

# Compute PSF values at each position
intensity_values = [psf(xi, yi) for yi in y, xi in x]

# Can be visualized with any plotting library, e.g.
# using CairoMakie
# heatmap(x, y, intensity_values, colormap=:viridis)
```

## Relationship to Airy Function

The Gaussian2D model can approximate the Airy disk pattern using the empirical relationship:

```math
\sigma \approx 0.42 \frac{\lambda}{\text{NA}}
```

where λ is the wavelength and NA is the numerical aperture.

## Performance Considerations

The Gaussian2D PSF is the most computationally efficient model available in MicroscopePSFs.jl, making it ideal for:

- Large-scale simulations
- Rapid prototyping
- Applications requiring real-time performance
- SMLM fitting algorithms where speed is critical

## Limitations

While computationally efficient, the Gaussian2D model has several limitations:

1. It doesn't account for diffraction effects that cause rings in the actual PSF
2. It can't model defocus or other aberrations
3. It's less accurate at larger distances from the center
4. It doesn't capture the true 3D nature of microscope PSFs

For applications requiring higher physical accuracy, consider using `Airy2D`, `Scalar3D`, or `Vector3D` models instead.