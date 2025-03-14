# Conventions

This page documents the conventions used throughout MicroscopePSFs.jl.

## Coordinate System

MicroscopePSFs.jl uses the following coordinate system conventions:

- **Physical units**: All distances are in micrometers (μm)
- **Origin**: The origin (0,0,0) is at the PSF center
- **Axes**:
  - x-axis: Lateral direction (positive right)
  - y-axis: Lateral direction (positive down)
  - z-axis: Axial direction (positive away from objective)

### Image Coordinates

When generating PSF images:

- Images are returned as Julia matrices with indices `[row, column]`
- This corresponds to `[y, x]` in the coordinate system
- The first element `[1, 1]` is the top-left pixel

## Normalization

PSF models use the following normalization conventions:

- **Intensity**: PSF intensity is normalized to integrate to 1.0 over the entire domain
- **Amplitude**: Complex amplitudes are normalized such that `|amplitude|²` gives the normalized intensity
- **Generated images**: PSF images are normalized to sum to 1.0 by default, unless specified otherwise

## Parameters

Common parameters across PSF models:

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in micrometers (μm)
- `n`: Refractive index of the sample medium
- `n_i`: Refractive index of the immersion medium

## Aberrations

Aberrations are represented using Zernike polynomials:

- Zernike coefficients use the Noll indexing scheme by default
- Coefficients are specified in wavelength units
- Positive z (defocus) corresponds to increasing the optical path length
- Note that this is consistent with a sample moving away from the objective

## Units

| Parameter   | Unit                |
|:------------|:--------------------|
| Wavelength  | Micrometers (μm)    |
| Position    | Micrometers (μm)    |
| Pixel size  | Micrometers (μm)    |
| Sigma       | Micrometers (μm)    |
| NA          | Dimensionless       |
| Refractive index | Dimensionless  |
| Zernike coefficients | Wavelengths |

## Performance Considerations

- Array dimensions are typically ordered as `[y, x]` for 2D and `[y, x, z]` for 3D
- For performance-critical code, provide pre-allocated arrays whenever possible
- PSF models are designed to be thread-safe for parallel evaluation

## Type Hierarchy

The package uses the following type hierarchy structure:

```
AbstractPSF
├── Gaussian2D
├── Airy2D
├── Scalar3DPSF
├── Vector3DPSF
└── SplinePSF
```

The `Gaussian2D` and `Airy2D` types are 2D PSF models, while `Scalar3DPSF` and `Vector3DPSF` are 3D PSF models. Each implements the appropriate methods for PSF evaluation in their respective dimensions.