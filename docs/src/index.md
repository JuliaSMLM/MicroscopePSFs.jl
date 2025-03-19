# MicroscopePSFs.jl

A Julia package for simulating microscope point spread functions (PSFs).

## Overview

MicroscopePSFs.jl provides a flexible and performant framework for calculating microscope point spread functions (PSFs) using various models, from simple 2D Gaussians to full 3D vectorial models. The package is designed for use in single-molecule localization microscopy (SMLM) and other fields requiring accurate optical simulations.

## Features

- Multiple PSF models: GaussianPSF, AiryPSF, ScalarPSF, VectorPSF, and SplinePSF
- Common interface for all PSF types with function-call syntax `psf(x, y, z)`
- Realistic camera pixel integration with `integrate_pixels(psf, camera, emitter)`
- Support for optical aberrations via Zernike polynomials

## Installation

```julia
using Pkg
Pkg.add("MicroscopePSFs")
```

## Quick Start

```jldoctest; output = false
using MicroscopePSFs

# Create a GaussianPSF
psf = GaussianPSF(0.15)  # σ = 150nm

# Evaluate at a specific position
intensity = psf(0.1, 0.2)  # Intensity at (x,y) = (0.1μm, 0.2μm)

# Generate a PSF image
x = range(-1, 1, length=101)  # μm
y = range(-1, 1, length=101)  # μm
img = [psf(xi, yi) for yi in y, xi in x]

# Create a simulated microscope image
nx = 20
ny = 20
pixel_size = 0.1 # micron
camera = IdealCamera(nx, ny, pixel_size)  # 20x20 pixels, 100nm size
emitter = Emitter2D(1.0, 1.0, 1000.0)               # At (1μm, 1μm) with 1000 photons
pixels = integrate_pixels(psf, camera, emitter)

# output
pixels
```

## Contents

See the navigation menu for detailed documentation on interfaces, conventions, and specific PSF types.