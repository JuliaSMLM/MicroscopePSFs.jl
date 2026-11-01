# VectorPSF

The `VectorPSF` model implements a comprehensive three-dimensional point spread function based on vectorial diffraction theory. It accounts for polarization effects, high-NA phenomena, dipole emission patterns, refractive index interfaces, and arbitrary aberrations.

## Mathematical Model

The VectorPSF is based on the vector Debye integral originally derived by Richards & Wolf and later extended for oriented dipole emission and stratified media by Sheppard–Török, Enderlein, and others. The expression for the electric field at the image plane is:

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

## Complex fields and the dipole frame

`dipole_fields(psf, x, y, z)` returns one `[Ex, Ey]` complex vector per mutually
incoherent dipole: one for a fixed orientation or three in x, y, z order for a
rotating dipole. Intensity is `sum(sum(abs2, E) for E in dipole_fields(...))`.
For a fixed dipole its only entry equals `amplitude`; rotating dipoles have
no single coherent amplitude.

Emitter z and dipole z share one frame: +z points from the coverslip into
the sample, away from the objective. `DipoleVector(1.0, 0, 1.0)` matches
lab dipole `[1, 0, 1]/√2`. Opposed objectives still mirror `pz`: construct the
bottom fixed dipole as
`DipoleVector(px, py, -pz)`. See [Multi-channel PSFs](multichannel.md) for
analyzers, the 4Pi mirror rule, and photon normalization.

The sign is anchored by Backlund, Lew, Backer, Sahl, and Moerner,
*ChemPhysChem* **15**, 587 (2014), [Figure 3D](https://pmc.ncbi.nlm.nih.gov/articles/PMC3992256/#F3),
reprinting Lew, Backlund, and Moerner, [*Nano Letters* **13**, 3967 (2013)](https://doi.org/10.1021/nl304359p).
Membrane-normal dipoles appear shifted laterally toward the cell edge on
both sides of focus, with larger shifts farther from focus. This statement
does not depend on a z-axis or defocus sign convention.

Versions 0.5.6 and earlier had the axial term with the wrong sign, so a fixed
dipole with both in-plane and axial components came out mirrored in `pz`.
Pure x, y, and z dipoles and the rotating dipole had correct intensities.
The `amplitude` and `dipole_fields` of a pure z dipole change sign in 0.5.7.
For a rotating dipole, `dipole_fields(psf, x, y, z)[3]` and
`pupil_at(psf, z)[3]` also change sign, with intensities unchanged.

A `VectorPSF` file stores filled pupils. A file below format version 1.1.0 has
no axial-sign marker. It was written by 0.5.6 or earlier, which used the old
axial dipole sign, or by an unreleased 0.5.7 development build. Some such
development builds wrote the corrected sign with format 1.0.0; which one wrote
the file cannot be told from the file. Mixed-dipole fields and coherent
interference can be affected. Loading preserves the stored fields exactly.
Rebuild it from its parameters to be sure.

Each `VectorPupilFunction` records its stored fields' `axial_sign_convention` as
`:legacy` or `:current`; the seven-argument constructor defaults to `:current`
and accepts that keyword. Unmarked pre-1.1.0 axial or rotating vector PSFs are
conservatively marked `:legacy` on load. Saving writes the convention per pupil,
so re-saving legacy fields in format 1.1.0 still warns on the next load, once per
file even inside wrappers and containers. `copy`, `pupil_at`, and `apply_defocus!`
preserve the convention. Refilling with `fill_vector_pupils!` or `update_pupils!`
recomputes the fields and sets `:current`. A `FourPiPSF` warns at construction
when its arms carry different conventions because the interference sign is wrong.
