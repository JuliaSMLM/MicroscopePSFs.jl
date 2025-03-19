# VectorPSF

The `VectorPSF` model implements a comprehensive three-dimensional point spread function based on vectorial diffraction theory. It accounts for polarization effects, high-NA phenomena, dipole emission patterns, refractive index interfaces, and arbitrary aberrations, providing the highest level of physical accuracy among the PSF models.

## Mathematical Model

The VectorPSF implements the full Richards-Wolf vector diffraction theory, accounting for the vectorial nature of light. The key formula for the electric field at the image plane is:

$$PSF(x - x_i, y - y_i, z_i, z_s) = \sum_{m=x,y} \sum_{n=p_x, p_y, p_z} |F[h(k_x, k_y)w_{m,n}A(k_x, k_y)e^{\iota 2\pi(k_x x_i+k_y y_i+k_{z\mathrm{med}} z_i-k_{z\mathrm{imm}}z_s)}]|^2$$

Where:
-  $(x_i, y_i, z_i)$ is the emitter position where $z_i$ represents the depth above the coverslip
-  $z_s$ is the distance the sample stage was moved from nominal focus
-  $F$ denotes the Fourier transform operation 
-  $h(k_x, k_y)$ is the complex pupil function incorporating aberrations 
-  $w_{m,n}$ represents electric field components at the pupil plane
-  $A(k_x, k_y)$ is the apodization factor for energy conservation 
-  $k_{z\mathrm{med}}$ and $k_{z\mathrm{imm}}$ are z-components of wave vectors in the sample and immersion media

## Constructor and Parameters

```julia
VectorPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
         base_pupil::Union{Nothing, PupilFunction}=nothing,
         base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
         n_medium::Real=1.33,
         n_coverslip::Real=1.52,
         n_immersion::Real=1.52,
         z_stage::Real=0.0, 
         grid_size::Integer=128)
```

### Required Parameters

- `nₐ`: Numerical aperture of the objective
- `λ`: Wavelength of light in microns
- `dipole`: Dipole orientation vector specified using `DipoleVector`

### Optional Parameters

- `base_pupil`: Optional base aberration pupil function (default: none)
- `base_zernike`: Zernike coefficients for aberrations (default: none)
- `n_medium`: Refractive index of the sample medium (default: 1.33, water)
- `n_immersion`: Refractive index of the immersion medium (default: 1.52, oil)
- `n_coverslip`: Refractive index of the coverslip (default: 1.52, glass)
- `z_stage`: Distance the sample stage was moved (μm) (default: 0.0)
- `grid_size`: Size of grid for pupil function (default: 128)

## Dipole Orientation Options

The VectorPSF supports both fixed and rotating dipole orientations:

### Fixed Dipole

```julia
# Create a VectorPSF with specific dipole orientation
dipole_xy = DipoleVector(0.707, 0.707, 0.0)  # 45° in XY plane
psf = VectorPSF(1.4, 0.532, dipole_xy, n_medium=1.33)
```

### Rotating Dipole (Isotropic Emission)

```julia
# No dipole specification required - uses incoherent average of x, y, z dipoles
psf = VectorPSF(1.4, 0.532, n_medium=1.33)
```

The rotating dipole model represents a freely rotating fluorophore by calculating the incoherent sum of three orthogonal dipole orientations.

## Key Features

- **Polarization Effects**: Models the vectorial nature of light propagation
- **Dipole Emission Patterns**: Accounts for the orientation of fluorescent dipoles
- **High-NA Accuracy**: Correctly handles the optical physics at high numerical apertures
- **Refractive Index Interfaces**: Models light propagation across media interfaces
- **Supercritical Angle Fluorescence**: Captures SAF effects for emitters near coverslip

## Refractive Index Interfaces

The VectorPSF accounts for refractive index mismatches:

```julia
# Water sample with oil immersion objective
psf = VectorPSF(
    1.4, 0.532, dipole_z,
    n_medium=1.33,   # Water sample
    n_immersion=1.52  # Oil immersion
)
```

## Supercritical Angle Fluorescence

A notable feature of the VectorPSF model is its ability to accurately capture Supercritical Angle Fluorescence (SAF), which occurs when light emitted by fluorophores near the coverslip is collected at angles exceeding the critical angle. SAF contribution decreases as the emitter moves away from the coverslip.

## Examples

```julia
# Create a basic VectorPSF with Z-oriented dipole
dipole_z = DipoleVector(0.0, 0.0, 1.0)
psf = VectorPSF(
    1.4,                # Numerical aperture
    0.532,              # Wavelength in microns
    dipole_z,           # Z-oriented dipole
    n_medium=1.33,      # Sample is water
    n_immersion=1.52   # Immersion oil
)

# Create a PSF with aberrations
zc = ZernikeCoefficients(15)
zc[11] = 0.5  # Add spherical aberration (normalized to rms = 1.0)
psf_aberrated = VectorPSF(
    1.4, 0.532, dipole_z,
    n_medium=1.33,
    base_zernike=zc
)
```

## Limitations

1. **Computational Intensity**: Most computationally demanding PSF model
2. **Parameter Complexity**: Requires understanding of multiple physical parameters
3. **Speed**: Significantly slower than other PSF models, especially with fine grid sizes

For standard usage patterns, camera integration, and comparison with other PSF types, see the [PSF Overview](overview.md).