# Spline PSF

The `SplinePSF` model provides an efficient representation of point spread functions using B-spline interpolation. It serves as a computational acceleration technique that can significantly speed up PSF evaluations for complex models like Vector3DPSF or Scalar3DPSF.

## Overview

B-spline interpolation allows MicroscopePSFs.jl to:

1. **Accelerate complex PSF calculations**: Pre-compute a complex PSF model once, then evaluate it quickly using interpolation
2. **Represent experimental PSFs**: Store measured PSFs from real microscopes
3. **Ensure high accuracy**: Cubic B-splines provide continuous second derivatives for accurate interpolation

## Use Cases

The most common usage pattern for SplinePSF is to accelerate computation with complex PSF models:

```julia
using MicroscopePSFs
using CairoMakie

# Create a computationally expensive PSF model (dipole along z-axis)
dipole_z = DipoleVector(0.0, 0.0, 1.0)
vector_psf = Vector3DPSF(
    1.4,                # Numerical aperture
    0.532,              # Wavelength
    dipole_z,           # Dipole orientation
    n_medium=1.33       # Sample medium refractive index
)

# Define coordinate sampling grid (higher density = better accuracy)
x_range = range(-2.0, 2.0, length=41)  # μm
y_range = range(-2.0, 2.0, length=41)  # μm
z_range = range(-2.0, 2.0, length=21)  # μm

# Create a spline representation (this step may take time)
spline_psf = SplinePSF(vector_psf, x_range, y_range, z_range)

# Now computations with the spline version are much faster
# Benchmark comparison
x, y, z = 0.5, -0.2, 0.3
@time vector_psf(x, y, z)    # Slow
@time spline_psf(x, y, z)    # Fast

# Visual comparison
fig = Figure(size=(900, 400))
z_positions = [-1.0, 0.0, 1.0]
for (i, z) in enumerate(z_positions)
    # Calculate PSF images at this z position
    intensity_vector = [vector_psf(xi, yi, z) for yi in y_range, xi in x_range]
    intensity_spline = [spline_psf(xi, yi, z) for yi in y_range, xi in x_range]
    
    # Plot orignal PSF
    ax1 = Axis(fig[1, i], aspect=DataAspect(), 
              title="Vector3DPSF (z=$(z)μm)",
              xlabel="x (μm)", ylabel="y (μm)")
    hm1 = heatmap!(ax1, x_range, y_range, intensity_vector, colormap=:viridis)
    
    # Plot spline version
    ax2 = Axis(fig[2, i], aspect=DataAspect(), 
              title="SplinePSF (z=$(z)μm)",
              xlabel="x (μm)", ylabel="y (μm)")
    hm2 = heatmap!(ax2, x_range, y_range, intensity_spline, colormap=:viridis)
end

fig
```

## Constructor Methods

The SplinePSF can be constructed in several ways:

### From Another PSF Model (Most Common)

```julia
# Create a spline version of any existing PSF
SplinePSF(psf::AbstractPSF, 
          x_range::AbstractRange,
          y_range::AbstractRange,
          z_range::AbstractRange;
          order::Integer=3)

# For 2D PSFs (Gaussian2D, Airy2D)
SplinePSF(psf::AbstractPSF, 
          x_range::AbstractRange,
          y_range::AbstractRange;
          order::Integer=3)
```

### From PSF Data (For Experimental PSFs)

```julia
# From a 3D PSF stack with coordinate ranges
SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
          x_range::AbstractRange,
          y_range::AbstractRange,
          z_range::AbstractRange;
          order::Integer=3)

# From a 2D PSF image
SplinePSF(psf_image::AbstractArray{<:Real,2}, 
          x_range::AbstractRange,
          y_range::AbstractRange;
          order::Integer=3)
```

## Examples

### Accelerating a Scalar3DPSF PSF

```julia
using MicroscopePSFs
using CairoMakie
using BenchmarkTools

# Create a Scalar3DPSF PSF with aberrations
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomials
add_spherical!(zc, 0.5)       # Add spherical aberration
add_astigmatism!(zc, 0.3)     # Add astigmatism

scalar_psf = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)

# Create a spline version for faster evaluation
x_range = range(-2.0, 2.0, length=81)
y_range = range(-2.0, 2.0, length=81)
z_range = range(-2.0, 2.0, length=41)

# This step takes time but only needs to be done once
spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)

# Compare speed for single point evaluation
@btime scalar_psf(0.5, 0.5, 0.1)
@btime spline_psf(0.5, 0.5, 0.1)

# Compare speed for camera integration
pixel_edges_x = collect(0:0.1:2.0)
pixel_edges_y = collect(0:0.1:2.0)
camera = IdealCamera(pixel_edges_x, pixel_edges_y)  # 20x20 camera
emitter = Emitter3D(1.0, 1.0, 0.0, 1000.0)          # Emitter with 1000 photons

@btime integrate_pixels($scalar_psf, $camera, $emitter)
@btime integrate_pixels($spline_psf, $camera, $emitter)

# Visualization of integrated pixels
pixels_scalar = integrate_pixels(scalar_psf, camera, emitter)
pixels_spline = integrate_pixels(spline_psf, camera, emitter)

# Camera physical coordinates
x_phys = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) / 2
y_phys = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) / 2

# Plot comparison
fig = Figure(size=(800, 400))

ax1 = Axis(fig[1, 1], aspect=DataAspect(), 
          title="Scalar3DPSF",
          xlabel="x (μm)", ylabel="y (μm)")
ax1.yreversed = true
hm1 = heatmap!(ax1, x_phys, y_phys, pixels_scalar', colormap=:viridis)
scatter!(ax1, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3], aspect=DataAspect(), 
          title="SplinePSF",
          xlabel="x (μm)", ylabel="y (μm)")
ax2.yreversed = true
hm2 = heatmap!(ax2, x_phys, y_phys, pixels_spline', colormap=:viridis)
scatter!(ax2, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
Colorbar(fig[1, 4], hm2)

fig
```

### Using Experimental PSF Data

```julia
using MicroscopePSFs
using CairoMakie
using HDF5

# Load experimental PSF data (example)
# In practice, you would load your own measured bead stack
function load_experimental_psf(filename)
    h5open(filename, "r") do file
        return read(file, "psf_stack")
    end
end

# Define physical coordinates for the experimental PSF
pixel_size = 0.1  # μm
z_step = 0.2      # μm
measured_psf = load_experimental_psf("path/to/measured_psf.h5")

# Create coordinate ranges based on PSF stack dimensions
nx, ny, nz = size(measured_psf)
x_range = range(-nx//2 * pixel_size, (nx//2-1) * pixel_size, length=nx)
y_range = range(-ny//2 * pixel_size, (ny//2-1) * pixel_size, length=ny)
z_range = range(-nz//2 * z_step, (nz//2-1) * z_step, length=nz)

# Create the SplinePSF
experimental_psf = SplinePSF(measured_psf, x_range, y_range, z_range)

# Use it like any other PSF
pixel_edges_x = collect(0:0.1:2.0)
pixel_edges_y = collect(0:0.1:2.0)
camera = IdealCamera(pixel_edges_x, pixel_edges_y)
emitter = Emitter3D(1.0, 1.0, 0.0, 1000.0)
pixels = integrate_pixels(experimental_psf, camera, emitter)
```

## B-Spline Interpolation Details

The SplinePSF uses cubic B-spline interpolation by default (order=3), which provides:

1. Continuous second derivatives for smooth interpolation
2. Local support for efficient computation
3. Good approximation of smooth functions
4. Analytic derivatives for optimization
5. Minimal memory requirements

## When to Use SplinePSF

Consider using SplinePSF when:

1. **Performance is critical**: For fitting or simulation tasks requiring many PSF evaluations
2. **Working with complex models**: To accelerate Vector3DPSF or aberrated Scalar3DPSF evaluations
3. **Using experimental PSFs**: When working with measured PSF data
4. **Running on limited compute resources**: To make complex models feasible on modest hardware