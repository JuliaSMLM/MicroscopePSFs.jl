# Airy2D

The `Airy2D` PSF model represents the diffraction-limited point spread function for a circular aperture under the paraxial approximation. It is more physically accurate than the Gaussian approximation, particularly for modeling the characteristic diffraction rings.

## Mathematical Model

The Airy pattern is defined in terms of its field amplitude:

```math
A(r) = \frac{\nu}{\sqrt{4\pi}} \cdot \frac{2J_1(\nu r)}{\nu r}
```

where:
- ``r = \sqrt{x^2 + y^2}`` is the radial distance from the optical axis
- ``\nu = \frac{2\pi \cdot \text{NA}}{\lambda}`` is the optical parameter
- ``J_1`` is the Bessel function of the first kind, order 1
- ``\text{NA}`` is the numerical aperture
- ``\lambda`` is the wavelength in microns

The intensity is given by:

```math
I(r) = |A(r)|^2
```

## Constructor

```julia
Airy2D(na, wavelength)
```

### Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns

### Examples

```julia
# Create an Airy PSF for a high-NA objective with green light
psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm

# Create from a Gaussian2D PSF (for comparison purposes)
gaussian_psf = Gaussian2D(0.15)
airy_equivalent = Airy2D(gaussian_psf, λ=0.532)
```

## Methods

### Evaluation

```julia
# Evaluate PSF at a specific position
intensity = psf(x, y)

# Get complex amplitude
amp = amplitude(psf, x, y)
```

### Creating Images

```julia
# Create a grid of positions
x_coords = -2:0.1:2  # microns
y_coords = -2:0.1:2  # microns

# Compute PSF values at each position
intensity_values = [psf(xi, yi) for yi in y_coords, xi in x_coords]

# Can be visualized with any plotting library, e.g.
# using CairoMakie
# heatmap(x_coords, y_coords, intensity_values, colormap=:viridis)
```

## Properties of the Airy Pattern

The Airy pattern has several notable features:

1. **Central Maximum**: Contains 83.8% of the total intensity
2. **First Minimum**: Occurs at a radius of 1.22λ/NA (Rayleigh criterion)
3. **First Ring**: Contains 7.2% of the total intensity
4. **Subsequent Rings**: Containing decreasing fractions of the intensity

## Performance Considerations

The Airy2D PSF is moderately computationally efficient:
- More expensive than Gaussian2D due to Bessel function evaluations
- Much faster than 3D models like Scalar3D and Vector3D
- Good balance between physical accuracy and performance

## Relationship to Other PSF Models

- **Gaussian Approximation**: The Airy pattern can be approximated by a Gaussian with σ ≈ 0.42λ/NA
- **3D Models**: The Airy2D pattern is the in-focus (z=0) slice of more complex 3D models

## Limitations

The Airy2D model has several limitations:

1. Only valid for in-focus imaging (no defocus modeling)
2. Uses paraxial approximation, less accurate for high-NA objectives
3. Doesn't account for aberrations, polarization effects, or refractive index mismatches
4. Limited to circular, uniform apertures

For applications requiring 3D imaging or higher physical accuracy, consider using `Scalar3D` or `Vector3D` models.