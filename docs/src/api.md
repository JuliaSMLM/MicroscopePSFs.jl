# API Reference

This page provides a comprehensive reference of the types and functions exported by MicroscopePSFs.jl.

## PSF Types

### 2D PSF Models

* `Gaussian2D` - Gaussian approximation PSF
* `Airy2D` - Airy disk PSF for diffraction-limited systems

### 3D PSF Models

* `Scalar3DPSF` - 3D PSF using scalar diffraction theory
* `Vector3DPSF` - 3D PSF using vectorial diffraction theory

### Data-driven PSF

* `SplinePSF` - B-spline representation of arbitrary PSFs

## Core Interface Functions

* `amplitude` - Calculate complex field amplitude
* `integrate_pixels` - Integrate PSF over camera pixels
* `integrate_pixels_amplitude` - Integrate complex amplitude over pixels

## Pupil Functions

* `PupilFunction` - Complex pupil function representation
* `VectorPupilFunction` - Vector pupil for polarization effects

## Zernike Aberrations

* `ZernikeCoefficients` - Representation of optical aberrations
* `add_defocus!` - Add defocus aberration
* `add_astigmatism!` - Add astigmatism aberration
* `add_coma!` - Add coma aberration
* `add_spherical!` - Add spherical aberration

## Emitters

* `DipoleVector` - 3D vector representing dipole orientation
* `DipoleEmitter3D` - 3D emitter with dipole orientation

Note: To create dipole vectors with specific orientations, use the `DipoleVector` constructor:

```julia
# X-oriented dipole
dipole_x = DipoleVector(1.0, 0.0, 0.0)

# Y-oriented dipole
dipole_y = DipoleVector(0.0, 1.0, 0.0)

# Z-oriented dipole
dipole_z = DipoleVector(0.0, 0.0, 1.0)
```

## I/O Functions

* `save_psf` - Save PSF to a file
* `load_psf` - Load PSF from a file
* `save_spline_psf` - Save SplinePSF to a file
* `load_spline_psf` - Load SplinePSF from a file