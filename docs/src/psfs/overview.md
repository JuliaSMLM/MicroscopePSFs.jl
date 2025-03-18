# PSF Types Overview

MicroscopePSFs.jl provides several PSF (Point Spread Function) models with varying complexity, accuracy, and computational requirements. This page gives an overview of the available models and guidance on selecting the appropriate PSF for your application.

## PSF Model Comparison

| PSF Type | Parameters | 2D/3D | Aberrations | Polarization | Relative Speed |
|:---------|:------------|:-------|:------------|:--------------|:---------------|
| `Gaussian2D` | σ | 2D | No | No | Fastest |
| `Airy2D` | NA, λ | 2D | No | No | Fast |
| `Scalar3DPSF` | NA, λ, n | 3D | Yes | No | Moderate |
| `Vector3DPSF` | NA, λ, dipole, n_medium, etc. | 3D | Yes | Yes | Slowest |
| `SplinePSF` | any | 2D/3D | Via source PSF | Via source PSF | Fast evaluation |

## When to Use Each PSF Type

### Gaussian2D

Best for:
- Rapid prototyping and initial development
- Simple fitting algorithms where computational speed is critical
- Applications where physical accuracy is less important than performance
- Educational purposes demonstrating basic PSF concepts

### Airy2D

Best for:
- Diffraction-limited 2D imaging simulations
- More accurate 2D fitting that accounts for diffraction rings
- Cases where you need a physically accurate model but don't need 3D capabilities
- Applications requiring a good balance between accuracy and speed

### Scalar3DPSF

Best for:
- 3D imaging simulations with moderate accuracy requirements
- Modeling defocus, spherical aberration, and other aberrations
- Applications with moderate NA objectives (typically NA < 1.2)
- When you need 3D capabilities but polarization effects aren't critical

### Vector3DPSF

Best for:
- High-NA objectives (NA > 1.2)
- Applications where polarization effects matter
- Modeling complex dipole emission patterns
- Research requiring the highest physical accuracy
- Simulations with significant refractive index mismatches
- When you need to account for all types of optical aberrations

### SplinePSF

Best for:
- Accelerating computationally intensive PSF models
- Using experimental PSF measurements from calibration beads
- Performance-critical applications like real-time fitting
- Creating fast approximations of complex physical models

## Standard Usage Pattern

All PSF types follow the same core interface, making it easy to switch between models:

```julia
# Create a PSF (example with Airy2D)
psf = Airy2D(1.4, 0.532)  # NA=1.4, λ=532nm

# Evaluate at specific position
intensity = psf(0.1, 0.2)  # at x=0.1μm, y=0.2μm

# Get complex field amplitude
amp = amplitude(psf, 0.1, 0.2)

# Create image grid
x = y = range(-1, 1, length=101)  # μm
img = [psf(xi, yi) for yi in y, xi in x]

# Camera integration example
pixel_size = 0.1  # μm
camera = IdealCamera(1:20, 1:20, pixel_size)
emitter = Emitter2D(1.0, 1.0, 1000.0)  # x, y, photons
pixels = integrate_pixels(psf, camera, emitter)
```

## Visual Comparison

Below is a comparison of the different PSF models using the same physical parameters:

```jldoctest; output = false 
using MicroscopePSFs
using CairoMakie

function compare_psf_profiles()
    # Common parameters
    na = 1.4
    wavelength = 0.532  # μm
    n = 1.518
    
    # Create consistent PSF instances
    gaussian = Gaussian2D(0.22 * wavelength / na)
    airy = Airy2D(na, wavelength)
    scalar = Scalar3DPSF(na, wavelength, n)
    
    # Dipole for vector PSF (z-oriented)
    dipole_z = DipoleVector(1.0, 0.0, 0.0)
    vector = Vector3DPSF(na, wavelength, dipole_z, n_medium=n)
    
    # Define positions for profile
    x = range(-1, 1, length=200)  # μm
    y = 0.0
    z = 0.0
    
    # Calculate profiles
    gaussian_profile = [gaussian(xi, y) for xi in x]
    airy_profile = [airy(xi, y) for xi in x]
    scalar_profile = [scalar(xi, y, z) for xi in x]
    vector_profile = [vector(xi, y, z) for xi in x]
    
    # Create figure
    fig = Figure(size=(900, 700))
    
    # 1D profiles
    ax1 = Axis(fig[1, 1:2], 
              xlabel="Position (μm)", 
              ylabel="Normalized Intensity",
              title="PSF Intensity Profiles (y=0)")
    
    lines!(ax1, x, gaussian_profile, label="Gaussian2D", linewidth=2, color=:blue)
    lines!(ax1, x, airy_profile, label="Airy2D", linewidth=2, color=:red)
    lines!(ax1, x, scalar_profile, label="Scalar3DPSF", linewidth=2, color=:green)
    lines!(ax1, x, vector_profile, label="Vector3DPSF", linewidth=2, color=:purple)
    
    axislegend(ax1, position=:rt)
    
    # 2D images
    x_grid = y_grid = range(-1, 1, length=101)  # μm
    
    # Calculate 2D images
    gaussian_img = [gaussian(xi, yi) for yi in y_grid, xi in x_grid]
    airy_img = [airy(xi, yi) for yi in y_grid, xi in x_grid]
    scalar_img = [scalar(xi, yi, 0.0) for yi in y_grid, xi in x_grid]
    vector_img = [vector(xi, yi, 0.0) for yi in y_grid, xi in x_grid]
    
    # Plot 2D images
    psf_images = [gaussian_img, airy_img, scalar_img, vector_img]
    titles = ["Gaussian2D", "Airy2D", "Scalar3DPSF", "Vector3DPSF"]
    
    for (i, (img, title)) in enumerate(zip(psf_images, titles))
        ax = Axis(fig[2, i], aspect=DataAspect(), 
                 title=title,
                 xlabel="x (μm)", 
                 ylabel=i==1 ? "y (μm)" : "")
        
        hm = heatmap!(ax, x_grid, y_grid, img, colormap=:viridis)
    end
    
    return fig
end

fig = compare_psf_profiles()

# output
Figure()
```

## PSF Conversion

The PSF types provide methods to convert between models where appropriate:

```julia
# Convert Airy2D to Gaussian2D approximation
airy = Airy2D(1.4, 0.532)
gaussian = Gaussian2D(airy)  # Automatically uses appropriate σ

# Convert Gaussian2D to equivalent Airy2D
gaussian = Gaussian2D(0.15)
airy = Airy2D(gaussian, λ=0.532)  # Need to specify wavelength
```

## Computational Considerations

When selecting a PSF model, consider the trade-off between accuracy and computational cost:

- **Gaussian2D**: Closed-form expression, extremely fast
- **Airy2D**: Uses Bessel functions, fast but more expensive than Gaussian
- **Scalar3DPSF**: Requires numerical integration or pre-computation, moderate speed
- **Vector3DPSF**: Most computationally intensive, especially with aberrations
- **SplinePSF**: Fast evaluation but requires initial computation or measurement

For performance-critical applications, consider using `SplinePSF` to pre-compute a more complex PSF model for faster repeated evaluations.