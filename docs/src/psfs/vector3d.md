# Vector3D

The `Vector3DPSF` model implements a comprehensive three-dimensional point spread function based on vectorial diffraction theory. This model accounts for polarization effects, high-NA phenomena, dipole emission patterns, and arbitrary aberrations, providing the highest level of physical accuracy.

## Mathematical Model

The Vector3D PSF uses the Richards-Wolf vector diffraction theory to calculate the electromagnetic field near focus:

```math
\mathbf{E}(\mathbf{r}) = \int_{pupil} \mathbf{P}(\theta, \phi) e^{i k (\mathbf{s} \cdot \mathbf{r})} \sin\theta \, d\theta \, d\phi
```

where:
- ``\mathbf{E}(\mathbf{r})`` is the electric field vector at position ``\mathbf{r} = (x, y, z)``
- ``\mathbf{P}(\theta, \phi)`` is the complex vector pupil function
- ``\mathbf{s}`` is the unit vector in the direction of wave propagation
- ``\theta`` and ``\phi`` are the polar and azimuthal angles in the pupil
- ``k = 2\pi / \lambda`` is the wave number
- The intensity is calculated as ``I(\mathbf{r}) = |\mathbf{E}(\mathbf{r})|^2``

The model accounts for polarization, refractive index mismatches, and emitter dipole orientation.

## Constructor

```julia
Vector3DPSF(na, wavelength, dipole; n_medium=1.33, n_immersion=1.518, n_coverslip=1.518,
          focal_z=0.0, base_zernike=nothing, zernike_coeffs=nothing, grid_size=128)
```

### Required Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns
- `dipole`: Dipole orientation vector specified using `DipoleVector`

### Optional Parameters

- `n_medium`: Refractive index of the sample medium (defaults to 1.33)
- `n_immersion`: Refractive index of the immersion medium (defaults to 1.518)
- `n_coverslip`: Refractive index of the coverslip (defaults to 1.518)
- `focal_z`: Focal plane position in microns (defaults to 0.0)
- `base_zernike`: Pre-created Zernike coefficients for aberrations (defaults to nothing)
- `zernike_coeffs`: Alternative name for Zernike coefficients (defaults to nothing)
- `grid_size`: Size of grid for pupil function (defaults to 128)

### Examples

```julia
# Create a basic Vector3D PSF with Z-oriented dipole
dipole_z = DipoleVector(0.0, 0.0, 1.0)
psf = Vector3DPSF(
    1.4,                # Numerical aperture
    0.532,              # Wavelength in microns
    dipole_z,           # Z-oriented dipole
    n_medium=1.33,      # Sample is water
    n_immersion=1.518   # Immersion oil
)

# Create a PSF with X-oriented dipole
dipole_x = DipoleVector(1.0, 0.0, 0.0)
psf_x = Vector3DPSF(
    1.4,
    0.532,
    dipole_x,
    n_medium=1.33
)

# Create a PSF with aberrations
zc = ZernikeCoefficients(15)
add_spherical!(zc, 0.5)
psf_aberrated = Vector3DPSF(
    1.4,
    0.532,
    dipole_z,
    n_medium=1.33,
    base_zernike=zc
)
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
x_coords = -2:0.1:2  # microns
y_coords = -2:0.1:2  # microns
z = 0.0  # focal plane

# Compute PSF values at each position
intensity_2d = [psf(xi, yi, z) for yi in y_coords, xi in x_coords]

# Create a 3D PSF stack
z_stack = range(-2, 2, length=5)  # z positions in microns
intensity_3d = Array{Float64}(undef, length(y_coords), length(x_coords), length(z_stack))
for (k, zi) in enumerate(z_stack)
    for (j, yi) in enumerate(y_coords)
        for (i, xi) in enumerate(x_coords)
            intensity_3d[j, i, k] = psf(xi, yi, zi)
        end
    end
end
```

## Dipole Orientation

The Vector3D PSF can model different emitter dipole orientations:

```julia
# Create dipole vectors for specific orientations
dipole_x = DipoleVector(1.0, 0.0, 0.0)  # X-oriented dipole
dipole_y = DipoleVector(0.0, 1.0, 0.0)  # Y-oriented dipole
dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Z-oriented dipole

# Custom dipole orientation
dipole_xy = DipoleVector(0.707, 0.707, 0.0)  # 45° in the XY plane
psf = Vector3DPSF(1.4, 0.532, dipole_xy, n_medium=1.33)
```

## Aberration Modeling

The Vector3D PSF supports the same aberration modeling as Scalar3D:

```julia
zc = ZernikeCoefficients(21)  # Support up to 21st Zernike mode

# Add common aberrations
add_defocus!(zc, 1.0)
add_astigmatism!(zc, 0.5)
add_coma!(zc, 0.3)
add_spherical!(zc, 0.2)
add_aberration!(zc, 9, 0.1)  # Add trefoil using OSA index 9

# Create PSF with these aberrations
psf = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33, base_zernike=zc)
```

## Performance Considerations

The Vector3DPSF is the most computationally intensive model:

- Significantly slower than other PSF models
- Performance scales with the resolution of the pupil function
- Can benefit from pre-computation strategies
- Recommended for applications where high physical accuracy is essential

## When to Use Vector3D

The Vector3DPSF model is recommended for:

1. High-NA objectives (NA > 1.2)
2. Polarization-sensitive applications
3. Accurate modeling of dipole emission patterns
4. Applications with significant refractive index mismatches
5. Research requiring the highest physical accuracy