# Conventions

This page documents the conventions used throughout MicroscopePSFs.jl.

## Coordinate System

MicroscopePSFs.jl uses the following coordinate system conventions:

- **Physical units**: All distances are in micrometers (μm)
- **Origin**: For PSF evaluation, the origin (0,0,0) is at the PSF center
- **Axes**:
  - x-axis: Lateral direction (positive right)
  - y-axis: Lateral direction (positive down)
  - z-axis: Axial direction (positive away from objective)

## Image and Camera Coordinates

When generating PSF images with `integrate_pixels`:

- Camera coordinates have (0,0) at the top-left corner of the top-left pixel
- Images are returned as Julia matrices with indices `[row, column]`
- This corresponds to `[y, x]` in the coordinate system
- The first element `[1, 1]` is the top-left pixel

## Normalization

PSF models use the following normalization conventions:

- **PSF evaluation**: PSF intensity is normalized to integrate to 1.0 over the entire $x,y$ plane
- **Amplitude**: Complex amplitudes are normalized such that `|amplitude|²` gives the normalized intensity
- **Camera images**: Images generated with `integrate_pixels` contain actual photon counts based on emitter.photons

## Aberrations

Aberrations are represented using Zernike polynomials:

- Zernike polynomials are normalized to give an RMS of 1.0
- Coefficients can be set directly via indexing: `zc[5] = 0.5`
- The OSA indexing scheme is used by default
- Positive z (defocus) corresponds to the emitter moving away from the objective

## Units

| Parameter   | Unit                |
|:------------|:--------------------|
| Wavelength  | Micrometers (μm)    |
| Position    | Micrometers (μm)    |
| Pixel size  | Micrometers (μm)    |
| Sigma       | Micrometers (μm)    |
| NA          | Dimensionless       |
| Refractive index | Dimensionless  |
| Zernike coefficients | Normalized (RMS=1) |

## Type Hierarchy

The package uses the following type hierarchy structure:

```
AbstractPSF
├── Abstract2DPSF
│   ├── GaussianPSF
│   └── AiryPSF
├── Abstract3DPSF
│   ├── ScalarPSF
│   └── VectorPSF
└── SplinePSF
```

## Performance Considerations

- For multi-emitter simulations, use the `support` parameter to limit computation to relevant regions
- Use `SplinePSF` to accelerate complex PSF models
- The package automatically uses multithreading in the `integrate_pixels` methods. 
