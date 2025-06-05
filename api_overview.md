# MicroscopePSFs.jl API Overview

This document provides a concise overview of the MicroscopePSFs.jl package API, designed for both human users and AI assistants.

## Key Concepts

MicroscopePSFs.jl models **point spread functions (PSFs)** for microscopy applications, particularly single-molecule localization microscopy (SMLM). The package provides:

- **2D and 3D PSF models** from simple Gaussian to complex vectorial diffraction
- **Aberration modeling** using Zernike polynomials
- **Dipole orientation support** for modeling fluorophore emission patterns
- **Camera integration** for simulating realistic microscope images
- **Performance optimization** through spline interpolation and multi-threading

**Units Convention**: All distances are in **microns** (μm). Wavelengths should be specified in microns (e.g., 0.532 for 532nm green light).

## Type Hierarchy

```
AbstractPSF{T}
├── Abstract2DPSF{T}
│   ├── GaussianPSF{T}         # Isotropic 2D Gaussian
│   └── AiryPSF{T}             # 2D Airy diffraction pattern
├── Abstract3DPSF{T}
│   ├── ScalarPSF{T}           # 3D scalar diffraction
│   └── VectorPSF{T}           # 3D vectorial with dipole
└── SplinePSF{T,IT}            # Spline interpolation (2D/3D)

PupilFunction{T}               # Complex pupil for scalar models
VectorPupilFunction{T}         # Vectorial pupil with polarization
```

## Essential Types

### 2D PSF Models

```julia
# Gaussian PSF - simplest model
struct GaussianPSF{T} <: Abstract2DPSF{T}
    σ::T        # Standard deviation in microns
end

# Airy PSF - diffraction-limited 2D
struct AiryPSF{T} <: Abstract2DPSF{T}
    nₐ::T       # Numerical aperture
    λ::T        # Wavelength in microns
end
```

### 3D PSF Models

```julia
# Scalar diffraction PSF - fast 3D model
struct ScalarPSF{T} <: Abstract3DPSF{T}
    nₐ::T                           # Numerical aperture
    λ::T                            # Wavelength in microns
    n::T                            # Refractive index of immersion medium
    pupil::PupilFunction{T}         # Complex pupil function
    # ... (internal fields)
end

# Vectorial PSF - accurate model with dipole orientation
struct VectorPSF{T} <: Abstract3DPSF{T}
    nₐ::T                           # Numerical aperture
    λ::T                            # Wavelength in microns
    dipole::DipoleVector{T}         # Dipole orientation
    pupil::VectorPupilFunction{T}   # Vector pupil function
    z_stage::T                      # Stage position offset
    # ... (internal fields)
end
```

### Support Types

```julia
# Dipole orientation (normalized vector)
struct DipoleVector{T}
    px::T
    py::T
    pz::T
end

# Zernike aberration coefficients
struct ZernikeCoefficients{T}
    phase::Vector{T}      # Phase coefficients (radians)
    mag::Vector{T}        # Magnitude coefficients
end
```

## Constructor Examples

### Basic PSF Construction

```julia
using MicroscopePSFs

# 2D Gaussian PSF with σ = 150nm
psf_gauss = GaussianPSF(0.15)

# 2D Airy PSF: NA=1.4, λ=532nm
psf_airy = AiryPSF(1.4, 0.532)

# 3D Scalar PSF: NA=1.4, λ=680nm, n=1.52 (oil immersion)
psf_scalar = ScalarPSF(1.4, 0.680, 1.52)

# Vector PSF with x-oriented dipole
dipole = DipoleVector(1.0, 0.0, 0.0)  # x-dipole
psf_vector = VectorPSF(1.4, 0.690, dipole; z_stage=0.0)
```

### PSF with Aberrations

```julia
# Create Zernike coefficients (max radial order 15)
zc = ZernikeCoefficients(15)

# Add specific aberrations (Noll indexing)
zc.phase[5] = 0.3    # Defocus
zc.phase[6] = 0.5    # Vertical astigmatism
zc.phase[7] = -0.2   # Oblique astigmatism

# Create aberrated PSF
psf_aberrated = ScalarPSF(1.4, 0.532, 1.52; zernike_coeffs=zc)
```

### Spline Interpolation

```julia
# Convert any PSF to spline for fast evaluation
x_range = -2.0:0.05:2.0  # microns
y_range = -2.0:0.05:2.0
z_range = -1.0:0.05:1.0

psf_spline = SplinePSF(psf_scalar, x_range, y_range, z_range)
```

## Core Functions

### PSF Evaluation

```julia
# 2D PSF evaluation
intensity = psf(x, y)                    # Single point
intensity = psf(x_vec, y_vec)            # Vectorized

# 3D PSF evaluation  
intensity = psf(x, y, z)                 # Single point
intensity = psf(x_vec, y_vec, z_vec)     # Vectorized

# Complex amplitude (before squaring)
amp = amplitude(psf, x, y)              # 2D
amp = amplitude(psf, x, y, z)           # 3D
```

### Image Generation

```julia
# Setup camera
camera = IdealCamera(nx, ny, pixel_size)  # nx, ny pixels, size in microns

# Create emitter
emitter = Emitter2D(x, y, photons)       # 2D
emitter = Emitter3D(x, y, z, photons)    # 3D

# Generate image
image = integrate_pixels(psf, camera, emitter)

# With support region (faster computation)
image = integrate_pixels(psf, camera, emitter; support=0.5)

# Multiple emitters
emitters = [Emitter3D(x[i], y[i], z[i], photons[i]) for i in 1:n]
image = integrate_pixels(psf, camera, emitters)
```

### Amplitude Integration

```julia
# Integrate complex amplitude (preserves phase)
complex_image = integrate_pixels_amplitude(psf, camera, emitter)

# Useful for coherent imaging or interferometry
```

## Common Workflows

### 1. Simple 2D Imaging

```julia
using MicroscopePSFs, SMLMData

# Create PSF
psf = GaussianPSF(0.15)  # 150nm standard deviation

# Setup camera: 100x100 pixels, 100nm pixel size
camera = IdealCamera(100, 100, 0.1)

# Create emitter at (5μm, 5μm) with 1000 photons
emitter = Emitter2D(5.0, 5.0, 1000.0)

# Generate image
image = integrate_pixels(psf, camera, emitter)
```

### 2. 3D Localization with Astigmatism

```julia
# Create PSF with astigmatism for 3D
zc = ZernikeCoefficients(5)
zc.phase[6] = 1.0  # Add astigmatism

psf = ScalarPSF(1.4, 0.680, 1.52; zernike_coeffs=zc)

# Emitter at different z positions
for z in [-0.5, 0.0, 0.5]  # microns
    emitter = Emitter3D(5.0, 5.0, z, 1000.0)
    image = integrate_pixels(psf, camera, emitter)
    # PSF shape changes with z due to astigmatism
end
```

### 3. Dipole Orientation Imaging

```julia
# Create dipoles with different orientations
dipole_x = DipoleVector(1.0, 0.0, 0.0)
dipole_z = DipoleVector(0.0, 0.0, 1.0)

psf_x = VectorPSF(1.4, 0.690, dipole_x)
psf_z = VectorPSF(1.4, 0.690, dipole_z)

# Compare PSF patterns
emitter = DipoleEmitter3D(5.0, 5.0, 0.0, 1000.0, dipole_x)
image_x = integrate_pixels(psf_x, camera, emitter)

emitter_z = DipoleEmitter3D(5.0, 5.0, 0.0, 1000.0, dipole_z)
image_z = integrate_pixels(psf_z, camera, emitter)
```

### 4. Performance Optimization with Splines

```julia
# Original PSF (slow for many evaluations)
psf = ScalarPSF(1.4, 0.680, 1.52)

# Convert to spline (fast evaluation)
psf_spline = SplinePSF(psf, -2:0.05:2, -2:0.05:2, -1:0.05:1)

# Use spline for faster multi-emitter simulations
emitters = [Emitter3D(rand()*10, rand()*10, randn()*0.5, 1000.0) 
            for _ in 1:100]
            
@time image_slow = integrate_pixels(psf, camera, emitters)
@time image_fast = integrate_pixels(psf_spline, camera, emitters)
# Spline version is typically 10-100x faster
```

## Complete Examples

### Example 1: Simulating STORM/PALM Data

```julia
using MicroscopePSFs, SMLMData, Random

# Microscope parameters
NA = 1.4                    # Oil objective
wavelength = 0.647          # 647nm in microns
n_immersion = 1.515         # Oil refractive index

# Create PSF
psf = ScalarPSF(NA, wavelength, n_immersion)

# Camera setup: 256x256 pixels, 160nm pixel size
camera = IdealCamera(256, 256, 0.160)

# Generate random emitters
Random.seed!(42)
n_emitters = 50
emitters = [
    Emitter3D(
        rand() * 40.0,      # x in [0, 40] microns
        rand() * 40.0,      # y in [0, 40] microns  
        randn() * 0.3,      # z in [-0.9, 0.9] microns
        500 + rand() * 2000 # 500-2500 photons
    )
    for _ in 1:n_emitters
]

# Generate image
image = integrate_pixels(psf, camera, emitters; support=0.7)

# Add noise (Poisson + Gaussian read noise)
# Note: Use your preferred noise simulation method
# For example: noisy_image = image + sqrt.(image) .* randn(size(image))  # Simple Poisson approximation
noisy_image = image .+ 0.01 * randn(size(image))  # Gaussian read noise
```

### Example 2: PSF Engineering for 3D Imaging

```julia
using MicroscopePSFs

# Create tetrapod PSF using phase mask
zc = ZernikeCoefficients(15)

# Tetrapod phase mask coefficients (example values)
zc.phase[11] = 1.5   # Primary spherical
zc.phase[22] = 0.8   # Secondary spherical
zc.phase[5] = 0.3    # Defocus

# Create engineered PSF
psf_tetrapod = ScalarPSF(1.4, 0.680, 1.52; zernike_coeffs=zc)

# Evaluate PSF at different z positions
z_positions = -1.5:0.1:1.5
camera = IdealCamera(64, 64, 0.1)

for (i, z) in enumerate(z_positions)
    emitter = Emitter3D(3.2, 3.2, z, 1000.0)
    image = integrate_pixels(psf_tetrapod, camera, emitter)
    # PSF rotates with z, enabling 3D localization
end
```

### Example 3: Analyzing PSF Properties

```julia
using MicroscopePSFs

# Create PSF
psf = ScalarPSF(1.4, 0.532, 1.52)

# Compute FWHM at focus
x = -1.0:0.01:1.0
y = 0.0
z = 0.0
intensity = [psf(xi, y, z) for xi in x]
intensity_norm = intensity ./ maximum(intensity)

# Find FWHM
half_max_indices = findall(intensity_norm .>= 0.5)
fwhm = x[half_max_indices[end]] - x[half_max_indices[1]]
println("Lateral FWHM: $(fwhm * 1000) nm")

# Axial profile
z_range = -2.0:0.01:2.0
axial_intensity = [psf(0.0, 0.0, zi) for zi in z_range]

# 3D OTF (via FFT of PSF)
# ... (implementation depends on specific analysis needs)
```

## Important Notes

1. **Memory efficiency**: Use `support` parameter in `integrate_pixels` to limit computation region
2. **Threading**: Multi-emitter integration is automatically parallelized
3. **Automatic differentiation**: Compatible with Enzyme.jl for gradient computations
4. **Normalization**: PSFs are normalized so total integrated intensity equals photon count
5. **Coordinate system**: (0,0) is at the corner of the first pixel, not the center

## Common Pitfalls

- Remember all distances are in **microns**, not nanometers
- Zernike coefficients use **Noll indexing** (starting from 1)
- For `VectorPSF`, the dipole vector is automatically normalized
- Support regions should be chosen carefully - too small truncates the PSF, too large wastes computation
- When using splines, ensure the interpolation range covers all possible emitter positions