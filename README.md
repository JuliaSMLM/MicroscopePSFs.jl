# MicroscopePSFs.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/MicroscopePSFs.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/MicroscopePSFs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/MicroscopePSFs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/MicroscopePSFs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/MicroscopePSFs.jl)

A Julia package for working with microscope Point Spread Functions (PSFs). This package provides implementations of common PSF models and tools for integrating them with camera geometry for single-molecule localization microscopy applications.

## Features

- Multiple PSF implementations (2D and 3D models)
- Integration with camera geometry via the SMLMData.jl package
- Complex field amplitude calculations for coherent optics
- Flexible pixel integration with adjustable sampling density
- Zernike polynomial tools for wavefront modeling

## PSF Models Overview

| PSF Type | Description | Key Parameters | When to Use |
|:---------|:------------|:---------------|:------------|
| `GaussianPSF` | Isotropic 2D Gaussian | `σ` | Rapid prototyping, computational efficiency |
| `AiryPSF` | Diffraction-limited circular aperture | `NA`, `λ` | Accurate 2D diffraction modeling |
| `ScalarPSF` | 3D scalar diffraction model | `NA`, `λ`, `n` | 3D imaging, aberrations |
| `VectorPSF` | 3D vectorial model with polarization | `NA`, `λ`, dipole, refractive indices | High-NA objectives, polarization effects |
| `SplinePSF` | B-spline approximation | Sampled grid | Accelerating computation of complex PSFs |

For detailed descriptions of each PSF type, see the [documentation](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable/psfs/overview/).

## Core Interface

All PSF types implement a consistent interface:

```julia
# Direct evaluation
intensity = psf(x, y)       # 2D position (x,y in μm)
intensity = psf(x, y, z)    # 3D position (3D PSFs only)

# Complex amplitude
amp = amplitude(psf, x, y)  # Returns complex field amplitude

# Camera integration
pixels = integrate_pixels(psf, camera, emitter)  # Realistic image formation
```

See the [interface documentation](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable/interface/) for details.

## Installation

```julia
using Pkg
Pkg.add("MicroscopePSFs")
```

## Basic Usage

```julia
using MicroscopePSFs

# Create a PSF model
psf = AiryPSF(1.4, 0.532)  # NA = 1.4, λ = 532nm
# or
psf = GaussianPSF(0.15)    # σ = 150nm

# Direct PSF evaluation at a point
intensity = psf(0.5, 0.3)  # x = 0.5μm, y = 0.3μm

# Get complex field amplitude
amp = amplitude(psf, 0.5, 0.3)

# Create camera and emitter
camera = IdealCamera(20, 20, 0.1)  # 20×20 pixels, 100nm pixel size
emitter = Emitter2D(1.0, 1.0, 1000.0)  # x = 1μm, y = 1μm, 1000 photons

# Generate realistic microscope image
pixels = integrate_pixels(psf, camera, emitter)
```

## Simulating a 3D PSF with Aberrations

```julia
# Create a 3D PSF with aberrations
zernike = ZernikeCoefficients(15)        # Up to 15 Zernike terms
zernike.phase[6] = 0.5  # Add vertical astigmatism
psf_3d = ScalarPSF(1.4, 0.532, 1.518; zernike_coeffs=zernike)

# Create 3D emitter
emitter_3d = Emitter3D(1.0, 1.0, 0.5, 1000.0)  # x, y, z, photons

# Generate image
pixels_3d = integrate_pixels(psf_3d, camera, emitter_3d)
```

## Coordinate Systems and Units

- All physical dimensions are in micrometers (μm)
- Function arguments are (x, y, z)
- Array dimensions are [y, x, z]
- Camera coordinates have (0,0) at the top-left corner of the top-left pixel

For detailed conventions, see the [documentation](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable/conventions/).

## Documentation

Comprehensive documentation is available:
- [Stable Release](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable)
- [Development Version](https://JuliaSMLM.github.io/MicroscopePSFs.jl/dev)

## License

This project is licensed under the MIT License - see the LICENSE file for details.