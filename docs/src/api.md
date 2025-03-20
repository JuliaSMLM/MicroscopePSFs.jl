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

#### PSF Acceleration

```@docs
SplinePSF
```

### Core Interface Functions

```@docs
amplitude
integrate_pixels
integrate_pixels_amplitude
```

### Multi-Emitter Integration

The following functions support integrating PSFs for multiple emitters:

```@docs
integrate_pixels(psf::AbstractPSF, camera::AbstractCamera, emitters::Vector{<:AbstractEmitter}; support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf, sampling::Integer=2, threaded::Bool=true)
integrate_pixels_amplitude(psf::AbstractPSF, camera::AbstractCamera, emitters::Vector{<:AbstractEmitter}; support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf, sampling::Integer=2, threaded::Bool=true)
```

### Pupil Functions

```@docs
PupilFunction
VectorPupilFunction
```

### Zernike Module

```@docs
MicroscopePSFs.Zernike.ZernikeCoefficients
```

### Zernike Polynomial Functions

```@docs
MicroscopePSFs.Zernike.zernikepolynomial
MicroscopePSFs.Zernike.radialpolynomial
MicroscopePSFs.Zernike.max_radial_order
MicroscopePSFs.Zernike.evaluate_pupil
```

### Zernike Analysis Functions

```@docs
MicroscopePSFs.Zernike.rms
MicroscopePSFs.Zernike.significant_terms
```

### Index Conversion

```@docs
MicroscopePSFs.Zernike.nl2osa
MicroscopePSFs.Zernike.osa2nl
MicroscopePSFs.Zernike.nl2noll
MicroscopePSFs.Zernike.noll2nl
MicroscopePSFs.Zernike.osa2noll
MicroscopePSFs.Zernike.noll2osa
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