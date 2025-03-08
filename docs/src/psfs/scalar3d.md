# Scalar3D

The `Scalar3D` PSF model implements a three-dimensional point spread function based on scalar diffraction theory. This model accounts for defocus and aberrations using a complex pupil function approach, providing a good balance between physical accuracy and computational efficiency.

## Mathematical Model

The Scalar3D PSF uses the Fourier optics approach to calculate the field and intensity distribution:

```math
U(\mathbf{r}) = \int_{pupil} P(\boldsymbol{\rho}) e^{i k \boldsymbol{\rho} \cdot \mathbf{r}} d\boldsymbol{\rho}
```

where:
- ``U(\mathbf{r})`` is the complex field amplitude at position ``\mathbf{r} = (x, y, z)``
- ``P(\boldsymbol{\rho})`` is the complex pupil function at pupil coordinates ``\boldsymbol{\rho}``
- ``k = 2\pi / \lambda`` is the wave number
- The intensity is calculated as ``I(\mathbf{r}) = |U(\mathbf{r})|^2``

The pupil function can incorporate various aberrations, typically represented using Zernike polynomials.

## Constructor

```julia
Scalar3DPSF(na, wavelength, n; pupil=nothing, pupil_data=nothing, coeffs=nothing)
```

### Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns
- `n`: Refractive index of the medium

### Optional Parameters

- `pupil`: Pre-created `PupilFunction` instance
- `pupil_data`: Complex matrix to initialize the pupil function
- `coeffs`: `ZernikeCoefficients` instance for representing aberrations

### Examples

```julia
# Create an unaberrated 3D PSF
psf = Scalar3DPSF(1.4, 0.532, 1.518)

# Create a PSF with spherical aberration
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial
add_spherical!(zc, 0.5)       # Add 0.5 waves of spherical aberration
psf_aberrated = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)
```

## Methods

### Evaluation

```julia
# Evaluate PSF at a specific 3D position
intensity = psf(x, y, z)

# Get complex amplitude
amp = amplitude(psf, x, y, z)
```

### Creating Images

```julia
# Create a grid of positions for a 2D slice at a specific z
x = range(-1, 1, length=31)
y = range(-1, 1, length=31)
z = 0.0 # focal plane

# Compute PSF values at each position
intensity_2d = [psf(xi, yi, z) for yi in y, xi in x]

# Create a 3D PSF stack
z_stack = range(-2, 2, length=5)  # z positions in microns
intensity_3d = Array{Float64}(undef, length(y), length(x), length(z_stack))
for (k, zi) in enumerate(z_stack)
    for (j, yi) in enumerate(y)
        for (i, xi) in enumerate(x)
            intensity_3d[j, i, k] = psf(xi, yi, zi)
        end
    end
end
```

## Aberration Modeling

The Scalar3D PSF can model optical aberrations through its pupil function:

```julia
# Create Zernike coefficients
zc = ZernikeCoefficients(15)  # Support up to 15th Zernike mode

# Add common aberrations
add_defocus!(zc, 1.0)         # 1 wave of defocus
add_astigmatism!(zc, 0.5)     # 0.5 waves of astigmatism
add_coma!(zc, 0.3)            # 0.3 waves of coma
add_spherical!(zc, 0.2)       # 0.2 waves of spherical aberration

# Create PSF with these aberrations
psf = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)
```

## Performance Considerations

The Scalar3D PSF offers several performance characteristics:

- More computationally intensive than 2D models (Gaussian2D, Airy2D)
- Significantly faster than the full vectorial model (Vector3D)
- Performance scales with the resolution of the pupil function
- Can be pre-computed for repeated evaluations

## Relationship to Other PSF Models

- **Airy2D**: The in-focus (z=0) slice of an unaberrated Scalar3D PSF is equivalent to the Airy pattern
- **Vector3D**: The Scalar3D is a simplified version of Vector3D that doesn't account for polarization effects

## Limitations

The Scalar3D model has several limitations:

1. Doesn't account for polarization effects, which become significant at high NA (>1.2)
2. May not accurately model certain high-NA phenomena
3. Simplified treatment of refractive index interfaces
4. Most accurate near the focal plane and for moderate NA values

For applications requiring higher physical accuracy with high-NA objectives, consider using the `Vector3D` model.