# Gaussian2D

The `Gaussian2D` PSF model represents the microscope point spread function as an isotropic 2D Gaussian function. While this is a mathematical approximation rather than a physical model derived from diffraction theory, it provides excellent computational efficiency for rapid prototyping and performance-critical algorithms.

## Mathematical Model

The Gaussian2D PSF is defined as:

```math
I(x, y) = \frac{1}{2\pi\sigma^2} \exp\left(-\frac{x^2 + y^2}{2\sigma^2}\right)
```

where:
- ``x, y`` are coordinates in physical units (microns)
- ``\sigma`` is the standard deviation in the same units

This function is normalized to integrate to 1 over the entire domain, ensuring energy conservation.

## Constructor and Parameters

```julia
Gaussian2D(σ::Real)
```

- `σ`: Standard deviation in microns, representing the width of the PSF

### Alternative Constructor
```julia
Gaussian2D(psf::Airy2D)  # Create from an Airy PSF
```

## Key Features

- **Computational Efficiency**: Fastest PSF model in the package, using a simple closed-form expression
- **Simplicity**: Simple mathematical form makes it ideal for prototyping
- **Approximation**: Provides a reasonable approximation of the central peak of diffraction-limited PSFs

## Examples

Creating a Gaussian PSF:

```julia
# Create a PSF with 150nm standard deviation
psf = Gaussian2D(0.15)

# Create a Gaussian approximation of an Airy disk
airy_psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm
gaussian_approximation = Gaussian2D(airy_psf)  # Automatically sets appropriate σ
```

## Relationship to Airy Function

The Gaussian2D model can approximate the Airy disk pattern using the empirical relationship:

```math
\sigma \approx 0.22 \frac{\lambda}{\text{NA}}
```

where λ is the wavelength and NA is the numerical aperture. This approximation works best near the center of the PSF.

## Limitations

- **No Diffraction Rings**: Doesn't capture the diffraction rings present in real microscope PSFs
- **No Defocus Modeling**: Can't model effects of defocus or 3D imaging
- **No Aberrations**: Doesn't account for optical aberrations
- **Simplified Physics**: Mathematical approximation rather than physically derived model
- **Less Accurate at Edges**: Diverges from physical PSFs at larger distances from the center

For standard usage patterns, camera integration, and comparison with other PSF types, see the [PSF Overview](overview.md).