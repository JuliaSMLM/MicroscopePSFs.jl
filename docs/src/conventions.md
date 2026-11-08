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
- Phase coefficients can be set with `zc[5, :phase] = 0.5` or `zc.phase[5] = 0.5`
- The Noll indexing scheme is used by default
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

## Numerical Range

The PSFs are checked over this range:

| Quantity | Range |
|:---------|:------|
| Wavelength | 0.1 to 10 μm |
| Numerical aperture | 0.01 to 1.7, and not above the refractive index |
| Refractive index | 1 to 4 |
| Gaussian σ | 0.01 to 10 μm |
| Positions and pixel edges | within 10⁴ μm (1 cm) of the origin |
| Pixel size | at least 10⁻³ μm (1 nm) |
| Axial position z | within 100 μm of focus |
| Photons per emitter | up to 10¹² |

Inside it, every value is finite and non-negative for `Float32`, `Float64` and integer
inputs and any mix of them. Integer pixel edges, positions and photon counts are converted
to floating point first, and give exactly the `Float64` result.

`Float32` keeps its own resolution. All-`Float32` images of an emitter near the origin agree
with `Float64` to within 0.1 percent of the image peak, and a `Float32` PSF evaluated at
`Float64` positions to within 0.04 percent. A `Float32` position 1 cm from the origin
resolves only about 1 nm, which moves the image of a PSF a few nanometers wide by up to a
few percent; use `Float64` positions there.

Outside the range, intermediate values such as squared wave numbers or pixel areas can
leave the range of the number type and give zero, infinity or `NaN`.

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
- The package automatically uses multithreading in the `integrate_pixels` methods by default
- For automatic differentiation, you can disable threading with the `threaded=false` parameter
- Disabling threading is required for compatibility with most automatic differentiation frameworks