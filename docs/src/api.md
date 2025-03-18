# API Reference

This page provides a comprehensive reference of the types and functions in MicroscopePSFs.jl.

## Public API

These types and functions form the core public API of MicroscopePSFs.jl.

### PSF Types

```@docs
AbstractPSF
```

#### 2D PSF Models

```@docs
GaussianPSF
AiryPSF
```

#### 3D PSF Models

```@docs
ScalarPSF
VectorPSF
```

#### Data-driven PSF

```@docs
SplinePSF
```

### Core Interface Functions

```@docs
amplitude
integrate_pixels
integrate_pixels_amplitude
```

### Pupil Functions

```@docs
PupilFunction
VectorPupilFunction
```

### Zernike Module

```@docs
ZernikeCoefficients
ZernikeIndexing
OSA
Noll
```

### Aberration Functions

```@docs
add_aberration!
add_defocus!
add_astigmatism!
add_coma!
add_spherical!
reset!
scale!
merge!
rms
trim!
significant_terms
```

### Zernike Polynomial Functions

```@docs
zernikepolynomial
radialpolynomial
max_radial_order
evaluate_pupil
```

### Index Conversion

```@docs
nl2osa
osa2nl
nl2noll
noll2nl
osa2noll
noll2osa
convert_index
```

### Emitters

```@docs
DipoleVector
DipoleEmitter3D
```

### I/O Functions

```@docs
save_psf
load_psf
```

## Complete API (All Documented Functions)

This section lists additional internal functions and types that are documented but not part of the public API.

```@autodocs
Modules = [MicroscopePSFs]
Public = false
```

### Zernike Module Internal API

```@autodocs
Modules = [MicroscopePSFs.Zernike]
Public = false
```