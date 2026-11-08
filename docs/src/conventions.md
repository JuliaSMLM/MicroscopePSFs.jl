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

Inside it, every PSF gives finite, non-negative intensities for `Float32`, `Float64` and
integer inputs and any mix of them.

Integer pixel edges, emitter positions (x, y and z), support radii and photon counts are
converted to floating point before any arithmetic. Fixed-width integers such as `Int`, `Int8`
and `UInt8` give exactly the result of the same numbers given as `Float64`; a `BigInt` gives
`BigFloat`. Integer optical parameters are accepted and follow Julia's promotion rules: with
integers and `Float64` the PSF is `Float64`, and an integer mixed with `Float32` gives a
`Float32` PSF.

A `SplinePSF` reproduces the samples it was built from. Between samples, cubic
interpolation can dip below zero, most where the samples are coarse compared with the width
of the PSF; the interpolated intensity is cut off at zero there.

`Float32` results carry `Float32` rounding, which grows with the size of the intermediate
values: the defocus phase, about 2π z n / λ, and positions far from the origin. For
example, an all-`Float32` scalar PSF at 100 μm of defocus differed from `Float64` by 0.24
percent of the image peak, and a `Float32` position 1 cm from the origin, which resolves
only about 1 nm, moved the image of a PSF a few nanometers wide by up to 2 percent. Use
`Float64` where that matters.

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