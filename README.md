# MicroscopePSFs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/MicroscopePSFs.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/MicroscopePSFs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/MicroscopePSFs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/MicroscopePSFs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/MicroscopePSFs.jl)

A Julia package for working with microscope Point Spread Functions (PSFs). This package provides implementations of common PSF models and tools for integrating them with camera geometry for single-molecule localization microscopy applications.

## Features

- Multiple PSF implementations (Gaussian2D, Airy2D, with Scalar3D and Vector3D coming soon)
- Integration with camera geometry via the SMLMData.jl package
- Complex field amplitude calculations for coherent optics
- Flexible pixel integration with adjustable sampling density
- Zernike polynomial tools for wavefront modeling

## Installation

```julia
using Pkg
Pkg.add("MicroscopePSFs")
```

## Basic Usage

```julia
using MicroscopePSFs

# Create a PSF model
psf = Airy2D(1.4, 0.532)  # NA = 1.4, λ = 532nm
# or
psf = Gaussian2D(0.15)    # σ = 150nm

# Direct PSF evaluation at a point
intensity = psf(0.5, 0.3)  # x = 0.5μm, y = 0.3μm

# Get complex field amplitude
amp = amplitude(psf, 0.5, 0.3)
```

## Coordinate Systems and Units

### Physical Units
All physical dimensions are in micrometers (μm):
- Positions (x, y, z)
- Wavelength (λ)
- PSF parameters (e.g., Gaussian σ)

### Coordinate Order
The package follows these coordinate ordering conventions:

1. Function arguments: Always (x, y, z)
```julia
intensity = psf(x, y)
amp = amplitude(psf, x, y)
```

2. Array dimensions: Always [y, x, z]
```julia
pixels = integrate_pixels(psf, camera, emitter)
# pixels has shape [ny, nx]
```

3. Field names: Use descriptive names (x, y, z)
```julia
emitter = Emitter2D(x=1.0, y=1.0, photons=1000)
```

### Camera Integration

The package integrates with camera geometry defined by the SMLMData.jl package. Pixel coordinates follow the standard image convention with (0,0) at the top-left:

```julia
# Create camera with 100 nm pixels
pixel_size = 0.1  # μm
nx, ny = 32, 32
x_edges = range(0, nx * pixel_size, nx + 1)
y_edges = range(0, ny * pixel_size, ny + 1)
camera = IdealCamera(x_edges, y_edges)

# Create emitter at (1μm, 1μm)
emitter = Emitter2D(1.0, 1.0, 1000.0)  # x, y, photons

# Integrate PSF over pixels
pixels = integrate_pixels(psf, camera, emitter, sampling=2)
# Returns [ny, nx] array of normalized intensities

# Get complex field amplitudes
amplitudes = integrate_pixels_amplitude(psf, camera, emitter, sampling=2)
# Returns [ny, nx] array of complex amplitudes
```

The `sampling` parameter controls subpixel sampling density for numerical integration. Higher values give more accurate results at the cost of computation time.

## PSF Models

### Gaussian2D
Represents an isotropic 2D Gaussian PSF:

```julia
psf = Gaussian2D(0.15)  # σ = 150nm
```

### Airy2D
Models a circular aperture under scalar, paraxial approximation:

```julia
psf = Airy2D(1.4, 0.532)  # NA = 1.4, λ = 532nm

# Convert between Gaussian and Airy
psf_gauss = Gaussian2D(psf)      # Approximate Airy as Gaussian
psf_airy = Airy2D(psf_gauss)    # Convert back (uses default λ)
psf_airy = Airy2D(psf_gauss, λ=0.488)  # Specify wavelength
```

### Scalar3D
Implements a scalar diffraction model for 3D PSFs:

```julia
# Create a 3D PSF with NA=1.4, λ=532nm, and refractive index n=1.518 (oil)
psf = Scalar3D(1.4, 0.532, 1.518)

# Evaluate PSF at a 3D position
intensity = psf(x, y, z)  # All coordinates in μm

# Get complex field amplitude in 3D
amp = amplitude(psf, x, y, z)

# Integrate over pixels at a specific z-plane
pixels = integrate_pixels(psf, camera, emitter3D)
```

The Scalar3D model accounts for defocus and spherical aberration effects, making it suitable for 3D imaging applications. It uses a physically accurate scalar diffraction model based on the Debye-Wolf formalism.

## Advanced Features

### Pixel Integration
The `integrate_pixels` function supports different sampling densities for accuracy vs. speed tradeoffs:

```julia
# Basic integration using pixel centers
pixels = integrate_pixels(psf, camera, emitter, sampling=1)

# More accurate integration with 4x4 subsampling
pixels = integrate_pixels(psf, camera, emitter, sampling=4)
```

### Complex Amplitudes
For coherent calculations, use the `amplitude` function and `integrate_pixels_amplitude`:

```julia
# Get complex field at a point
amp = amplitude(psf, x, y)

# Integrate complex amplitude over pixels
amps = integrate_pixels_amplitude(psf, camera, emitter)
intensities = abs2.(amps)  # Convert to intensity if needed
```

## Coming Soon

- Vector3D PSF model for high-NA systems with polarization effects
- GPU acceleration support for faster computation
- Additional aberration models and analysis tools

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests on our GitHub repository.

## License

This project is licensed under the MIT License - see the LICENSE file for details.