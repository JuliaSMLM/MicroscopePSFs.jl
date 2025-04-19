# MicroscopePSFs API Overview

This document provides a concise reference of the MicroscopePSFs Julia package API for use by LLMs when developing dependent packages.

## Core Concepts

MicroscopePSFs models point spread functions (PSFs) for microscopy applications, with:

- **Hierarchy of PSF models** (from simple to complex):
  - 2D Analytical (Gaussian, Airy)
  - 3D Scalar diffraction models
  - 3D Vector diffraction models with dipole emitters
  - Spline-based interpolated PSFs

- **Physical parameters** use standard microscopy conventions:
  - Wavelength (λ) in nanometers
  - Numerical aperture (NA) as dimensionless ratio
  - Refractive indices for media (typically oil, water, etc.)
  - Coordinates in same units (typically nanometers)

- **Emitter handling** supports both simple point emitters and oriented dipoles

- **Integration** capabilities for generating realistic camera images

## Exported Types

### PSF Types

```julia
AbstractPSF                      # Base abstract type for all PSFs

# 2D PSFs
GaussianPSF{T<:AbstractFloat}    # 2D Gaussian PSF with σ parameter 
AiryPSF{T<:AbstractFloat}        # 2D Airy PSF with NA and wavelength

# 3D PSFs
ScalarPSF{T<:AbstractFloat}      # 3D scalar diffraction theory PSF
VectorPSF{T<:AbstractFloat}      # 3D vectorial diffraction with dipoles
SplinePSF{T,IT}                  # Spline-interpolated PSF
```

### Pupil Functions

```julia
PupilFunction{T<:AbstractFloat}       # Complex-valued pupil function
VectorPupilFunction{T<:AbstractFloat} # Vector pupil function for polarization
```

### Aberration Modeling

```julia
ZernikeCoefficients               # Container for Zernike aberration coefficients
```

### Emitters (Re-exported from SMLMData)

```julia
AbstractEmitter                   # Base type for all point sources
Emitter2D                         # 2D point emitter (x,y)
Emitter3D                         # 3D point emitter (x,y,z)
DipoleEmitter3D                   # 3D dipole emitter with orientation
```

### Camera Models (Re-exported from SMLMData)

```julia
AbstractCamera                    # Base type for cameras
IdealCamera                       # Perfect photon detection
```

### Dipole Representation

```julia
DipoleVector                      # 3D dipole orientation
```

## Key Functions

### PSF Evaluation

```julia
# Evaluate PSF at coordinates
(psf::AbstractPSF)(x, y)                # 2D evaluation at position (x,y)
(psf::AbstractPSF)(x, y, z)             # 3D evaluation at position (x,y,z)

# Get complex field amplitude
amplitude(psf::AbstractPSF, x, y)       # 2D complex amplitude
amplitude(psf::AbstractPSF, x, y, z)    # 3D complex amplitude
```

### PSF Integration

```julia
# Integrate PSF over camera pixels
integrate_pixels(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    support=2,                          # Support region size in PSF widths
    sampling=10,                        # Supersampling factor 
    threaded=true                       # Enable multithreading
)

# Lower-level integration with explicit pixel boundaries
integrate_pixels(
    psf::AbstractPSF,
    pixel_edges_x::Vector,
    pixel_edges_y::Vector,
    emitter::AbstractEmitter;
    support=2,
    sampling=10,
    threaded=true
)

# Integrate complex field amplitude
integrate_pixels_amplitude(
    psf::AbstractPSF,
    camera::AbstractCamera,
    emitter::AbstractEmitter;
    support=2,
    sampling=10,
    threaded=true
)
```

### Zernike Functions

```julia
# From the Zernike module
zernikepolynomial(n, m, ρ, θ)     # Evaluate Zernike polynomial
radialpolynomial(n, m, ρ)         # Evaluate radial part

# Index conversion between different schemes
nl2osa(n, l)                      # Convert (n,l) to OSA/ANSI index
osa2nl(j)                         # Convert OSA/ANSI index to (n,l)
nl2noll(n, l)                     # Convert (n,l) to Noll index
noll2nl(j)                        # Convert Noll index to (n,l)
osa2noll(j)                       # Convert OSA/ANSI index to Noll
noll2osa(j)                       # Convert Noll index to OSA/ANSI
```

### I/O Functions

```julia
save_psf(filename, object; metadata=Dict()) # Save PSF to HDF5 file
load_psf(filename)                          # Load PSF from HDF5 file
```

### Dipole-Related Functions

```julia
calculate_pupil_field(...)        # Calculate pupil field for dipole emitters
calculate_fresnel_coefficients(...) # Calculate Fresnel coefficients
```

## Common Usage Patterns

1. **Creating a PSF model:**
   ```julia
   # Simple 2D model
   psf = GaussianPSF(sigma=100.0)  # 100nm standard deviation

   # Airy model for diffraction-limited imaging
   psf = AiryPSF(;λ=600.0, NA=1.4)  # 600nm emission, 1.4 NA objective

   # 3D scalar model
   psf = ScalarPSF(;λ=600.0, NA=1.4, n=1.518, zrange=(-1000.0, 1000.0, 21))
   
   # 3D vector model with dipole orientation
   psf = VectorPSF(;λ=600.0, NA=1.4, n=1.518, zrange=(-1000.0, 1000.0, 21))
   ```

2. **Evaluating at a point:**
   ```julia
   # Evaluate intensity at a given position
   intensity = psf(x, y)       # 2D
   intensity = psf(x, y, z)    # 3D
   ```

3. **Generating a camera image:**
   ```julia
   # Setup camera and emitter
   camera = IdealCamera(pixel_size=100.0)  # 100nm pixels
   emitter = Emitter3D(x=1000.0, y=1000.0, z=0.0, I=1000.0)
   
   # Generate image
   img = integrate_pixels(psf, camera, emitter)
   ```

4. **Adding aberrations:**
   ```julia
   # Create aberration coefficients (OSA/ANSI indexing)
   coefs = zeros(15)
   coefs[4] = 0.5  # 0.5 waves of defocus
   coefs[8] = 0.2  # 0.2 waves of spherical aberration
   
   # Create PSF with aberrations
   psf = ScalarPSF(;λ=600.0, NA=1.4, n=1.518, zernike_coefficients=coefs)
   ```

5. **Saving/loading PSFs:**
   ```julia
   save_psf("my_psf.h5", psf, metadata=Dict("description" => "My custom PSF"))
   loaded_psf = load_psf("my_psf.h5")
   ```