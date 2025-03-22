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

Current I/O version: v1.0.0