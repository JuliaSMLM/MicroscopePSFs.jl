# MicroscopePSFs.jl

A Julia package for simulating microscope point spread functions (PSFs).

## Overview

MicroscopePSFs.jl provides a flexible and performant framework for calculating microscope point spread functions (PSFs) using various models, from simple 2D Gaussians to full 3D vectorial models. The package is designed for use in single-molecule localization microscopy (SMLM) and other fields requiring accurate optical simulations.

## Features

- Multiple PSF models: Gaussian2D, Airy2D, Scalar3DPSF, Vector3DPSF, and Spline PSFs
- Common interface for all PSF types with function-call syntax `psf(x, y, z)`
- Realistic camera pixel integration with `integrate_pixels(psf, camera, emitter)`
- Support for optical aberrations via Zernike polynomials
- GPU acceleration (for supported models)
- Optimized for performance in single-molecule fitting applications
- Interoperability with the broader Julia ecosystem

## Installation

```julia
using Pkg
Pkg.add("MicroscopePSFs")
```

## Quick Start

```julia
using MicroscopePSFs
using CairoMakie

# Create a Gaussian2D PSF
psf = Gaussian2D(0.15)  # σ = 150nm

# Evaluate at a specific position
intensity = psf(0.1, 0.2)  # Intensity at (x,y) = (0.1μm, 0.2μm)

# Generate a PSF image
x = range(-1, 1, length=101)  # μm
y = range(-1, 1, length=101)  # μm
img = [psf(xi, yi) for yi in y, xi in x]

# Visualize the PSF
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="Gaussian PSF (σ=150nm)",
          xlabel="x (μm)", 
          ylabel="y (μm)")
hm = heatmap!(ax, x, y, img, colormap=:viridis)
Colorbar(fig[1, 2], hm)

# Create a simulated microscope image
pixel_edges_x = collect(0:0.1:2.0)  # Convert to Vector
pixel_edges_y = collect(0:0.1:2.0)  # Convert to Vector
camera = IdealCamera(pixel_edges_x, pixel_edges_y)  # 20x20 pixels, 100nm size
emitter = Emitter2D(1.0, 1.0, 1000.0)               # At (1μm, 1μm) with 1000 photons
pixels = integrate_pixels(psf, camera, emitter)

# Camera physical coordinates
x_phys = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) / 2
y_phys = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) / 2

# Visualize the camera image
ax2 = Axis(fig[2, 1:2], aspect=DataAspect(),
           title="Integrated Camera Image",
           xlabel="x (μm)", 
           ylabel="y (μm)")
ax2.yreversed = true  # Flip y axis to match camera convention
hm2 = heatmap!(ax2, x_phys, y_phys, pixels', colormap=:viridis)
scatter!(ax2, [emitter.x], [emitter.y], 
         color=:red, marker=:cross, markersize=15)

fig
```

## Contents

See the navigation menu for detailed documentation on interfaces, conventions, and specific PSF types.