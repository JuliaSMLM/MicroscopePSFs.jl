# Gaussian2D

The `Gaussian2D` PSF model represents the microscope point spread function as an isotropic 2D Gaussian function. While this is a mathematical approximation rather than a physical model derived from diffraction theory, it works well for many applications and is computationally very efficient, making it an excellent choice for rapid prototyping and performance-critical algorithms.

## Mathematical Model

The Gaussian2D PSF is defined as:

```math
I(x, y) = \frac{1}{2\pi\sigma^2} \exp\left(-\frac{x^2 + y^2}{2\sigma^2}\right)
```

where:
- ``x, y`` are coordinates in physical units (microns)
- ``\sigma`` is the standard deviation in the same units

This function is normalized to integrate to 1 over the entire domain, which ensures energy conservation.

## Constructor Options

```julia
Gaussian2D(σ::Real)
```

### Parameters
- `σ`: Standard deviation in microns

### Type Parameters
- `T`: Numeric precision type, automatically determined from input

## Basic Usage

### Creating a PSF

```julia
# Create a Gaussian PSF with 150nm standard deviation
psf = Gaussian2D(0.15)

# Create a Gaussian approximation of an Airy disk
airy_psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm
gaussian_approximation = Gaussian2D(airy_psf)  # Automatically sets appropriate σ
```

### Direct Evaluation

```julia
# Evaluate PSF at a specific position
intensity = psf(0.1, 0.2)  # At position x=0.1μm, y=0.2μm

# Get complex amplitude (returns sqrt of intensity for Gaussian)
amp = amplitude(psf, 0.1, 0.2)
```

### Creating a PSF Image

```julia
# Create a grid of positions
x = range(-1, 1, length=101)  # μm
y = range(-1, 1, length=101)  # μm

# Compute PSF values on the grid
intensity_values = [psf(xi, yi) for yi in y, xi in x]

# Visualize with CairoMakie
using CairoMakie
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="Gaussian PSF (σ=150nm)",
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

## Relationship to Airy Function

The Gaussian2D model can approximate the Airy disk pattern using the empirical relationship:

```math
\sigma \approx 0.22 \frac{\lambda}{\text{NA}}
```

where λ is the wavelength and NA is the numerical aperture. This approximation works best near the center of the PSF and becomes less accurate at the edges where the Airy pattern has diffraction rings.

## Performance Considerations

The Gaussian2D PSF is the most computationally efficient model available in MicroscopePSFs.jl, making it ideal for:

- Large-scale simulations with many PSFs
- Real-time applications
- Initial algorithm development and testing
- SMLM fitting algorithms where speed is critical

The model uses a simple closed-form expression that can be evaluated very efficiently, with no need for numerical integration or special functions (unlike other PSF models).

## Limitations

While computationally efficient, the Gaussian2D model has several limitations:

1. **No Diffraction Rings**: It doesn't account for the diffraction rings that appear in real microscope PSFs
2. **No Defocus Modeling**: It can't model the effects of defocus or 3D imaging
3. **No Aberrations**: It doesn't account for optical aberrations
4. **Simplified Physics**: It's a mathematical approximation rather than a model derived from physics
5. **Less Accurate at Edges**: It diverges from physical PSFs at larger distances from the center

## Relationship to Other PSFs

The Gaussian2D model is related to other PSF models in the following ways:

1. **Airy2D**: The Gaussian2D can be created from an Airy2D using `Gaussian2D(airy_psf)` which will set σ appropriately
2. **Airy2D Conversion**: A Gaussian2D can be converted to an approximate Airy2D using `Airy2D(gaussian_psf, λ=wavelength)`
3. **3D Models**: The Gaussian2D can be viewed as an in-focus slice of a 3D Gaussian PSF model (though not directly provided by this package)

For applications requiring higher physical accuracy, consider using `Airy2D`, `Scalar3DPSF`, or `Vector3DPSF` models instead.