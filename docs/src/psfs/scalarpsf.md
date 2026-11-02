# ScalarPSF

The `ScalarPSF` model implements a three-dimensional point spread function based on scalar diffraction theory. This model accounts for defocus and optical aberrations using a complex pupil function approach, providing a good balance between physical accuracy and computational efficiency.

## Mathematical Model

The ScalarPSF uses the Fourier optics approach to calculate the complex field distribution:

```math
U(\mathbf{r}) = \int_{pupil} P(\boldsymbol{\rho}) e^{i k \boldsymbol{\rho} \cdot \mathbf{r}} d\boldsymbol{\rho}
```

where:
- ``U(\mathbf{r})`` is the complex field amplitude at position ``\mathbf{r} = (x, y, z)``
- ``P(\boldsymbol{\rho})`` is the complex pupil function at pupil coordinates ``\boldsymbol{\rho}``
- ``k = 2\pi / \lambda`` is the wave number
- The intensity is calculated as ``I(\mathbf{r}) = |U(\mathbf{r})|^2``

The pupil function can incorporate various aberrations, typically represented using Zernike polynomials.

## Constructor and Parameters

```julia
ScalarPSF(nₐ::Real, λ::Real, n::Real; 
         pupil::Union{Nothing, PupilFunction}=nothing,
         pupil_data::Union{Nothing, AbstractMatrix}=nothing,
         zernike_coeffs::Union{Nothing, ZernikeCoefficients}=nothing,
         grid_size::Integer=128)
```

### Required Parameters

- `nₐ`: Numerical aperture of the objective
- `λ`: Wavelength of light in microns
- `n`: Refractive index of the medium

### Optional Parameters

- `pupil`: Pre-created `PupilFunction` instance
- `pupil_data`: Complex matrix to initialize the pupil function
- `zernike_coeffs`: `ZernikeCoefficients` instance for representing aberrations
- `grid_size`: Generated pupil grid size (default: 128). A supplied `pupil` or
  `pupil_data` retains its own grid.

Supply one of `pupil`, `pupil_data`, or `zernike_coeffs`. If none is provided,
an unaberrated pupil is created using `grid_size`.

## Pupil sampling and field of view

The pupil field repeats in x and y with period
`L = (grid_size - 1) * λ / (2nₐ)` in microns, using the actual pupil grid size.
For rendering, keep `support` below `L/2` or use a frame narrower than `L`.
Increasing `grid_size` widens `L` and improves defocused intensity accuracy.
For NA 1.4, wavelength 0.6 microns, and grid 128, `L ≈ 27.21` microns.

For NA 1.4, wavelength 0.6 microns, index 1.518, and 0.5 rad RMS astigmatism
(`zc.phase[6] = 0.5`), the following maximum absolute intensity errors were
measured over a 31×31 grid on `[-1.5, 1.5]` microns in x and y. Errors are divided
by the peak intensity of a 512-sample pupil reference at each z. These measured
values depend on the optics, aberrations, and sampled positions.

| `grid_size` | z = 0.5 μm | z = 1 μm | z = 2 μm |
|---:|---:|---:|---:|
| 64 | 2.32% | 2.90% | 9.25% |
| 128 | 0.53% | 0.95% | 2.16% |
| 256 | 0.18% | 0.22% | 0.58% |

### Evaluation cost and a defocused-image example

The default grid of 128 takes about four times the evaluation time of grid 64.
With NA 1.2, wavelength 0.520 microns, index 1.3, and an unaberrated pupil, a
21×21-pixel frame on `[-1, 1]` microns per axis, `support=3`, `sampling=10`, and
10,000 photons gives the following maximum absolute pixel errors divided by
the peak of a 1024-sample pupil reference. These finite-reference measurements
use the same pixel quadrature at each grid size.

| `grid_size` | z = 0 μm | abs(z) = 0.5 μm | abs(z) = 1 μm |
|---:|---:|---:|---:|
| 64 | 0.6245% | 4.5750% | 2.2149% |
| 128 | 0.1682% | 1.2189% | 0.8065% |

Over z = -1, -0.5, 0, 0.5, and 1 microns, warmed single-image medians from three
repetitions per z are 2.449–2.613 s at grid 64 and 9.811–10.269 s at grid 128.
An 81×81×21 PSF tabulation with NA 1.2, wavelength 0.6 microns, index 1.33, and
`zc.phase[6]=0.5` rad takes 27.80 s at grid 64 and 93.05 s at grid 128, including
its first tabulation and compilation. These are local measurements on a shared
host with four Julia threads; they are not portable latency guarantees.

## Key Features

- **3D Imaging**: Models PSF behavior throughout 3D space, not just in focus
- **Aberration Modeling**: Supports arbitrary optical aberrations via Zernike polynomials
- **Complex Field**: Provides access to both amplitude and phase information
- **Physical Realism**: Based on physical principles of scalar diffraction theory

## Aberration Modeling

A key feature of the ScalarPSF is its ability to incorporate optical aberrations using Zernike polynomials with Noll indexing and RMS normalization:

```julia
# Create Zernike coefficients object
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial (Noll indexed)

# Add common aberrations by directly setting coefficients
# Note: The package uses Noll indexing and RMS normalization
zc.phase[4] = 0.5   # 0.5 rad RMS of defocus (Noll index 4)
zc.phase[5] = 0.3   # 0.3 rad RMS of astigmatism (Noll index 5)
zc.phase[7] = 0.2   # 0.2 rad RMS of coma (Noll index 7)
zc.phase[11] = 0.1  # 0.1 rad RMS of spherical aberration (Noll index 11)

# Create PSF with these aberrations
psf = ScalarPSF(1.4, 0.532, 1.518, zernike_coeffs=zc)
```

## Examples

Creating a ScalarPSF:

```julia
# Create a basic unaberrated 3D PSF
psf = ScalarPSF(1.4, 0.532, 1.518)  # NA=1.4, λ=532nm, n=1.518

# Create a PSF with spherical aberration
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial (Noll indexed)
zc.phase[11] = 0.5            # Add 0.5 rad RMS of spherical aberration (Noll index 11)
psf_aberrated = ScalarPSF(1.4, 0.532, 1.518, zernike_coeffs=zc)

# Create a PSF with a pre-computed pupil function
pupil = PupilFunction(1.4, 0.532, 1.518, zc)
psf_from_pupil = ScalarPSF(1.4, 0.532, 1.518, pupil=pupil)
```

## Limitations

1. **No Polarization Effects**: Doesn't account for polarization, which becomes significant at high NA (>1.2)
2. **Simplified Refractive Index Interfaces**: Simplified treatment of refractive index interfaces
3. **Scalar Approximation**: Uses scalar diffraction theory instead of full vector theory
4. **Moderate NA Assumption**: Most accurate for moderate NA values and regions near the focal plane

For standard usage patterns, camera integration, and comparison with other PSF types, see the [PSF Overview](overview.md).