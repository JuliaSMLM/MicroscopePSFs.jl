# PSF Types Overview

MicroscopePSFs.jl provides several PSF (Point Spread Function) models with varying complexity and accuracy. This page provides an overview of the available models and guidance on selecting the appropriate model for your application.

## Available PSF Types

| PSF Type | Dimensions | Physics | Best For | Performance |
|:---------|:-----------|:--------|:---------|:------------|
| `Gaussian2D` | 2D | Simple approximation | Simple modeling, fitting | Fastest |
| `Airy2D` | 2D | Diffraction-limited | Accurate in-focus PSF | Fast |
| `Scalar3D` | 3D | Scalar diffraction | 3D PSF modeling | Moderate |
| `Vector3D` | 3D | Full vectorial model | High-NA, polarization effects | Slowest |
| `SplinePSF` | Any | Data-driven | Fast evaluation of other PSFs | Fast evaluation |

## Choosing a PSF Model

### Gaussian2D

The simplest PSF model, approximating the microscope PSF as a 2D Gaussian function. This is a mathematical approximation rather than a physical model, but it's computationally efficient and often sufficient for simple applications.

**Best for**: 
- Rapid prototyping
- Simple SMLM fitting
- Cases where speed is prioritized over accuracy

### Airy2D

Models the in-focus PSF of a diffraction-limited system using the Airy disk formula. More accurate than a Gaussian for in-focus imaging, particularly at the edges of the PSF.

**Best for**:
- Diffraction-limited imaging simulations
- More accurate 2D SMLM fitting
- Educational purposes demonstrating diffraction

### Scalar3D

A 3D PSF model based on scalar diffraction theory, accounting for defocus and spherical aberrations. Provides a good balance between accuracy and computational efficiency.

**Best for**:
- 3D imaging simulations
- Modeling defocus and spherical aberrations
- Applications where vectorial effects are negligible

### Vector3D

The most comprehensive physical model, implementing full vectorial diffraction theory. Accounts for polarization effects, high-NA phenomena, and all types of aberrations.

**Best for**:
- High-NA objectives (>1.2)
- Polarization-sensitive applications
- Situations requiring high physical accuracy
- Modeling complex aberrations

### SplinePSF

A computational acceleration technique that represents any PSF model using B-spline interpolation for faster evaluation. Create an accurate but slow model once, then use the spline version for repeated calculations.

**Best for**:
- Accelerating computationally intensive PSF models
- Representing experimental PSF measurements
- Performance-critical applications

## Computational Considerations

When selecting a PSF model, consider the trade-off between accuracy and computational cost:

- **Gaussian2D**: Extremely fast, closed-form expression
- **Airy2D**: Fast, using Bessel functions
- **Scalar3D**: Moderate speed, requiring numerical integration or pre-computation
- **Vector3D**: Most computationally intensive, especially with aberrations
- **SplinePSF**: Fast evaluation after initial computation/measurement

For large-scale simulations or fitting applications, the simpler models may be preferable unless you specifically need the effects modeled by the more complex approaches.

## PSF Model Comparison

```julia
using MicroscopePSFs
using CairoMakie

# Define microscope parameters
na = 1.4
wavelength = 0.532  # μm
n = 1.518

# Create PSF models
gaussian = Gaussian2D(sigma=wavelength/(2*na))
airy = Airy2D(na=na, wavelength=wavelength)
scalar = Scalar3DPSF(na=na, wavelength=wavelength, n=n)
vector = Vector3DPSF(na=na, wavelength=wavelength, n=n)

# Create position arrays
x = range(-2, 2, length=200)  # μm
y = 0.0

# Generate PSF profiles
profiles = [
    [psf(xi, y) for xi in x] for psf in [gaussian, airy, scalar, vector]
]

# Create plot
fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], 
    xlabel="Position (μm)", 
    ylabel="Normalized Intensity",
    title="PSF Model Comparison")

lines = []
colors = [:blue, :red, :green, :purple]
labels = ["Gaussian2D", "Airy2D", "Scalar3D", "Vector3D"]

for (i, profile) in enumerate(profiles)
    lines!(ax, x, profile, color=colors[i], linewidth=2)
end

Legend(fig[1, 2], 
    [LineElement(color=c) for c in colors],
    labels,
    "PSF Models")

# Create PSF images for comparison
x_img = y_img = range(-1, 1, length=101)
imgs = [[psf(xi, yi) for yi in y_img, xi in x_img] for psf in [gaussian, airy, scalar, vector]]

# Display PSF images
for (i, img) in enumerate(imgs)
    ax_img = Axis(fig[2, i], aspect=DataAspect(),
                 title=labels[i],
                 xlabel="x (μm)",
                 ylabel="y (μm)")
    heatmap!(ax_img, x_img, y_img, img, colormap=:viridis)
end

fig
```