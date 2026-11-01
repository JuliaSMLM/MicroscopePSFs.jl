# I/O Functionality

MicroscopePSFs.jl provides functions for saving and loading PSF objects to and from HDF5 files. These functions allow you to:

- Save PSF models for later use
- Share PSF models with others
- Store experimental or simulated PSFs for reuse
- Save custom pupil functions or Zernike coefficients

## Core Functions

### save_psf

```julia
save_psf(filename::String, object; metadata::Dict=Dict())
```

Save a PSF or related object to an HDF5 file.

**Arguments**
- `filename`: Path where the PSF will be saved
- `object`: Object to save (PSF, ZernikeCoefficients, PupilFunction, etc.)
- `metadata`: Optional dictionary of additional metadata to include

**Returns**
- `filename` for chaining

### load_psf

```julia
load_psf(filename::String)
```

Load a PSF or related object from an HDF5 file.

**Arguments**
- `filename`: Path to the HDF5 file containing a saved PSF object

**Returns**
- The loaded object with its original type (PSF, ZernikeCoefficients, PupilFunction, etc.)

## Vector pupil convention

Starting with the 0.5.7 release, saved files have `io_version = "1.1.0"`.
This version identifies the corrected axial dipole sign; no per-group marker
is needed. The minimum supported file version remains 1.0.0.

Loading a file below version 1.1.0 with a nonzero fixed dipole z component or a
rotating dipole emits one warning for the entire file, including vectors nested
inside wrappers and containers. Files without a version are treated as 1.0.0.
The file has no axial-sign marker: it was written by 0.5.6 or earlier, which
used the old axial dipole sign, or by an unreleased 0.5.7 development build.
Some such development builds wrote the corrected sign with format 1.0.0;
which one wrote the file cannot be told from the file. Mixed-dipole fields
and coherent interference can be affected. Loading preserves the stored fields
exactly. Rebuild it from its optical, dipole, and aberration parameters to be sure.

## Basic Usage

```julia
using MicroscopePSFs

# Create a PSF model
psf = AiryPSF(1.4, 0.532)  # NA = 1.4, λ = 532nm

# Save PSF to an HDF5 file
save_psf("my_airy_psf.h5", psf)

# Load PSF from file
loaded_psf = load_psf("my_airy_psf.h5")
```

## Metadata and Additional Information

You can add metadata when saving a PSF:

```julia
# Add custom metadata
metadata = Dict(
    "created_by" => "John Doe",
    "description" => "Standard Airy PSF for green light",
    "experiment_id" => "E123"
)

save_psf("my_airy_psf_with_metadata.h5", psf; metadata=metadata)
```

## Compatible Types

The I/O system supports saving and loading the following types:

- All PSF models (`GaussianPSF`, `AiryPSF`, `ScalarPSF`, `VectorPSF`, `SplinePSF`)
- Pupil functions (`PupilFunction`, `VectorPupilFunction`)
- Zernike coefficients (`ZernikeCoefficients`)

## Examples

### Saving and Loading Different PSF Types

```julia
# Create a GaussianPSF
gaussian_psf = GaussianPSF(0.15)
save_psf("gaussian_psf.h5", gaussian_psf)

# Create a ScalarPSF with aberrations
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5  # Add vertical astigmatism
scalar_psf = ScalarPSF(1.4, 0.532, 1.518; zernike_coeffs=zc)
save_psf("scalar_psf_with_aberrations.h5", scalar_psf)

# Create a SplinePSF from a ScalarPSF and save it
spline_psf = SplinePSF(scalar_psf, lateral_range=1.0, lateral_step=0.05, axial_step=0.1)
save_psf("accelerated_spline_psf.h5", spline_psf)
```

### Saving and Loading Zernike Coefficients

```julia
# Create Zernike coefficients
zc = ZernikeCoefficients(15)
zc.phase[4] = 0.5  # Add defocus 
zc.phase[11] = 0.2  # Add spherical aberration

# Save coefficients
save_psf("zernike_coefficients.h5", zc)

# Load coefficients
loaded_zc = load_psf("zernike_coefficients.h5")

# Use loaded coefficients to create a PSF
psf = ScalarPSF(1.4, 0.532, 1.518; zernike_coeffs=loaded_zc)
```

### Loading in a Different Session

```julia
using MicroscopePSFs
using HDF5

# Load a previously saved PSF
psf = load_psf("scalar_psf_with_aberrations.h5")

# The loaded PSF behaves exactly like the original
pixels = integrate_pixels(psf, camera, emitter)
```

## File Format Details

PSFs are saved in HDF5 format with a standardized structure:

- `/attributes`: Metadata and type information
  - `io_version`: Version of the I/O system
  - `psf_type`: Type of the stored object
  - `creation_date`: When the file was created
  - Custom metadata provided by the user
- `/parameters`: Physical parameters of the PSF
- `/data`: PSF-specific data (fields, grid values, etc.)
- `/zernike_coefficients` (optional): Zernike coefficients for aberrations

This format is designed to be portable and future-proof.

## Version Compatibility

The I/O system uses versioning to ensure forward and backward compatibility. Files created with newer versions of MicroscopePSFs.jl may not be readable by older versions.

Current I/O version: v1.1.0

## Multi-channel and wrapper layout

`save_psf` and `load_psf` also support `PlanePSF`, `AnalyzerPSF`, `FourPiPSF`,
and `MultiChannelPSF`, at the root or recursively in child groups. Root metadata
keywords work as before. The existing types keep their dataset layout.
`IO_VERSION` is `1.1.0`. No PSFLearning or Reactant
import is needed to load and evaluate a container.

A container has root attribute `psf_type = "MultiChannelPSF"` and this layout
(the registration size is in Julia dimension order):

```text
parameters/
    n_channels
    reference
    transmission         # length N
    z_offset             # length N
    registration         # size (2, 3, N)
channels/
    channel_1/           # attribute psf_type, then the channel's own layout
    ...
    channel_N/
```

Each wrapper stores its scalar parameter and complete child PSFs:

| Type | Parameter dataset | Child groups |
|:-----|:------------------|:-------------|
| `PlanePSF` | `parameters/dz` | `psf` |
| `AnalyzerPSF` | `parameters/component` | `psf` |
| `FourPiPSF` | `parameters/phase` | `top`, `bottom` |

Every child group has a bare type-name `psf_type` attribute. Loading uses the
public constructors, enforcing the mirror, field, registration, and photon
fraction requirements. Homogeneous channel types are narrowed, so a saved
collection of matching splines reloads with a concrete spline element type.

```julia
using MicroscopePSFs

p = ScalarPSF(1.2, 0.68, 1.33)
mc = MultiChannelPSF(; channels=[p, PlanePSF(p; dz=0.2)])
mktempdir() do dir
    file = joinpath(dir, "channels.h5")
    save_psf(file, mc; metadata=Dict("description" => "Dual plane"))
    restored = load_psf(file)
    restored.channels[2](0.0, 0.0, 0.1)
end
```
