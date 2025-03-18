# AiryPSF

The `AiryPSF` model represents the diffraction-limited point spread function for a circular aperture under the paraxial approximation. Unlike the simpler Gaussian approximation, this model accurately captures the characteristic diffraction rings that appear in real microscope images.

## Mathematical Model

The Airy pattern is defined in terms of its field amplitude:

```math
A(r) = \frac{\nu}{\sqrt{4\pi}} \cdot \frac{2J_1(\nu r)}{\nu r}
```

where:
- ``r = \sqrt{x^2 + y^2}`` is the radial distance from the optical axis in microns
- ``\nu = \frac{2\pi \cdot \text{NA}}{\lambda}`` is the optical parameter
- ``J_1`` is the Bessel function of the first kind, order 1
- ``\text{NA}`` is the numerical aperture
- ``\lambda`` is the wavelength in microns

The intensity is given by the squared magnitude of the amplitude:

```math
I(r) = |A(r)|^2
```

## Constructor and Parameters

```julia
AiryPSF(na::Real, wavelength::Real)
```

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns

### Alternative Constructor

```julia
AiryPSF(psf::GaussianPSF; λ::Real=0.532)
```

Creates an `AiryPSF` that approximates the provided `GaussianPSF`, using the specified wavelength.

## Key Features

- **Physical Accuracy**: Correctly models the diffraction pattern from a circular aperture
- **Diffraction Rings**: Captures the characteristic rings in real microscope images
- **Rayleigh Criterion**: Naturally demonstrates the Rayleigh resolution criterion
- **Computational Efficiency**: More physically accurate than Gaussian while still being computationally efficient

## Properties of the Airy Pattern

The Airy pattern has several notable physical properties:

1. **Central Maximum**: Contains 83.8% of the total intensity
2. **First Minimum**: Occurs at a radius of 1.22λ/NA (the Rayleigh criterion)
3. **First Ring**: Contains 7.2% of the total intensity
4. **Subsequent Rings**: Contain decreasing fractions of the intensity

## Examples

Creating and using an Airy PSF:

```julia
# Create an Airy PSF for a high-NA objective with green light
psf = AiryPSF(1.4, 0.532)  # NA=1.4, wavelength=532nm

# Create from a GaussianPSF
gaussian_psf = GaussianPSF(0.15)
airy_equivalent = AiryPSF(gaussian_psf, λ=0.532)
```

## Limitations

1. **2D Only**: Only valid for in-focus imaging (no defocus modeling)
2. **Paraxial Approximation**: Less accurate for very high-NA objectives (> 1.4)
3. **No Aberrations**: Doesn't account for optical aberrations
4. **No Polarization**: Doesn't model polarization effects
5. **No Refractive Index Mismatches**: Assumes uniform media

For standard usage patterns, camera integration, and comparison with other PSF types, see the [PSF Overview](overview.md).