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
ScalarPSF(na::Real, wavelength::Real, n::Real; 
         pupil::Union{Nothing, PupilFunction}=nothing,
         pupil_data::Union{Nothing, AbstractMatrix}=nothing,
         coeffs::Union{Nothing, ZernikeCoefficients}=nothing)
```

### Required Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns
- `n`: Refractive index of the medium

### Optional Parameters

- `pupil`: Pre-created `PupilFunction` instance
- `pupil_data`: Complex matrix to initialize the pupil function
- `coeffs`: `ZernikeCoefficients` instance for representing aberrations

You should provide exactly one of the optional parameters. If none are provided, an unaberrated pupil is created.

## Key Features

- **3D Imaging**: Models PSF behavior throughout 3D space, not just in focus
- **Aberration Modeling**: Supports arbitrary optical aberrations via Zernike polynomials
- **Complex Field**: Provides access to both amplitude and phase information
- **Physical Realism**: Based on physical principles of scalar diffraction theory

## Aberration Modeling

A key feature of the ScalarPSF is its ability to incorporate optical aberrations:

```julia
# Create Zernike coefficients object
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial

# Add common aberrations
add_defocus!(zc, 1.0)         # 1 wave of defocus
add_astigmatism!(zc, 0.5)     # 0.5 waves of astigmatism
add_coma!(zc, 0.3)            # 0.3 waves of coma
add_spherical!(zc, 0.2)       # 0.2 waves of spherical aberration

# Create PSF with these aberrations
psf = ScalarPSF(1.4, 0.532, 1.518, coeffs=zc)
```

## Examples

Creating a Scalar3DPSF:

```julia
# Create a basic unaberrated 3D PSF
psf = ScalarPSF(1.4, 0.532, 1.518)  # NA=1.4, λ=532nm, n=1.518

# Create a PSF with spherical aberration
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial
add_spherical!(zc, 0.5)       # Add 0.5 waves of spherical aberration
psf_aberrated = ScalarPSF(1.4, 0.532, 1.518, coeffs=zc)

# Create a PSF with a pre-computed pupil function
pupil = PupilFunction(1.4, 0.532, 1.518, zernike_coeffs)
psf_from_pupil = ScalarPSF(1.4, 0.532, 1.518, pupil=pupil)
```

## Limitations

1. **No Polarization Effects**: Doesn't account for polarization, which becomes significant at high NA (>1.2)
2. **Simplified Refractive Index Interfaces**: Simplified treatment of refractive index interfaces
3. **Scalar Approximation**: Uses scalar diffraction theory instead of full vector theory
4. **Moderate NA Assumption**: Most accurate for moderate NA values and regions near the focal plane

For standard usage patterns, camera integration, and comparison with other PSF types, see the [PSF Overview](overview.md).