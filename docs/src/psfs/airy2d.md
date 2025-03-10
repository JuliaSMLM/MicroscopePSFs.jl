# Airy2D

The `Airy2D` PSF model represents the diffraction-limited point spread function for a circular aperture under the paraxial approximation. Unlike the simpler Gaussian approximation, this model accurately captures the characteristic diffraction rings that appear in real microscope images, making it more physically accurate while still maintaining good computational efficiency.

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

## Constructor Options

```julia
Airy2D(na::Real, wavelength::Real)
```

### Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns

### Alternative Constructor

```julia
Airy2D(psf::Gaussian2D; λ::Real=0.532)
```

Creates an `Airy2D` PSF that approximates the provided `Gaussian2D` PSF, using the specified wavelength.

### Type Parameters
- `T`: Numeric precision type, automatically determined from input

## Basic Usage

### Creating a PSF

```julia
# Create an Airy PSF for a high-NA objective with green light
psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm

# Create from a Gaussian2D PSF (for comparison purposes)
gaussian_psf = Gaussian2D(0.15)
airy_equivalent = Airy2D(gaussian_psf, λ=0.532)
```

### Direct Evaluation

```julia
# Evaluate PSF at a specific position
intensity = psf(0.1, 0.2)  # At position x=0.1μm, y=0.2μm

# Get complex amplitude
amp = amplitude(psf, 0.1, 0.2)
```

### Creating a PSF Image

```julia
# Create a grid of positions
x = range(-2, 2, length=201)  # μm
y = range(-2, 2, length=201)  # μm

# Compute PSF values on the grid
intensity_values = [psf(xi, yi) for yi in y, xi in x]

# Visualize with CairoMakie
using CairoMakie
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="Airy PSF (NA=1.4, λ=532nm)",
          xlabel="x (μm)", ylabel="y (μm)")
hm = heatmap!(ax, x, y, intensity_values, colormap=:viridis)
Colorbar(fig[1, 2], hm)
fig
```

## Integration with Camera

```julia
# Create camera with 100nm pixels (20×20 pixel grid)
pixel_size = 0.1  # μm
camera = IdealCamera(1:20, 1:20, pixel_size)

# Create emitter at position (1μm, 1μm) with 1000 photons
emitter = Emitter2D(1.0, 1.0, 1000.0)

# Integrate PSF over pixels with 2×2 subsampling
pixels = integrate_pixels(psf, camera, emitter, sampling=2)

# Visualize the camera image
using CairoMakie
fig = Figure(size=(500, 400))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="Integrated Camera Image",
          xlabel="x (μm)", ylabel="y (μm)")
ax.yreversed = true  # Flip y-axis to match camera convention

# Get physical coordinates of pixel centers
x_centers = (1:20) * pixel_size .- pixel_size/2
y_centers = (1:20) * pixel_size .- pixel_size/2

hm = heatmap!(ax, x_centers, y_centers, pixels', colormap=:viridis)
scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
Colorbar(fig[1, 2], hm)
fig
```

## Properties of the Airy Pattern

The Airy pattern has several notable features that distinguish it from the Gaussian approximation:

1. **Central Maximum**: Contains 83.8% of the total intensity
2. **First Minimum**: Occurs at a radius of 1.22λ/NA (the Rayleigh criterion)
3. **First Ring**: Contains 7.2% of the total intensity
4. **Subsequent Rings**: Contain decreasing fractions of the intensity
5. **Infinite Extent**: Unlike the Gaussian, which asymptotically approaches zero, the Airy pattern extends to infinity with alternating rings

## Performance Considerations

The Airy2D PSF offers a good balance between physical accuracy and computational efficiency:

- More computationally intensive than Gaussian2D due to Bessel function evaluations
- Much faster than 3D models like Scalar3DPSF and Vector3DPSF
- Handles the special case at r=0 efficiently by using a series expansion
- Well-suited for applications where diffraction rings are important but 3D effects are not

## Limitations

The Airy2D model has several limitations:

1. **2D Only**: Only valid for in-focus imaging (no defocus modeling)
2. **Paraxial Approximation**: Less accurate for very high-NA objectives (> 1.4)
3. **No Aberrations**: Doesn't account for optical aberrations
4. **No Polarization**: Doesn't model polarization effects
5. **No Refractive Index Mismatches**: Assumes uniform media

## Relationship to Other PSFs

The Airy2D pattern is related to other PSF models in the following ways:

1. **Gaussian2D**: The Airy pattern can be approximated by a Gaussian with σ ≈ 0.22λ/NA
2. **Scalar3DPSF**: The Airy2D pattern is the in-focus (z=0) slice of the Scalar3DPSF with no aberrations
3. **Vector3DPSF**: For low-NA objectives, the in-focus Vector3DPSF approaches the Airy pattern

For applications requiring 3D imaging or higher physical accuracy, consider using `Scalar3DPSF` or `Vector3DPSF` models instead.