# MicroscopePSFs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/MicroscopePSFs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/MicroscopePSFs.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/MicroscopePSFs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/MicroscopePSFs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/MicroscopePSFs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/MicroscopePSFs.jl)

*MicroscopePSFs* provides a Microscope Point Spread Function (PSF) calculator.  

Current implementaions provide widefield PSFs assuming an incoherent point source.  

## Features

- 2D Gaussian PSF
- 2D Airy PSF
- 2D/3D Scalar PSF
- Scalar PSF allows arbitrary Pupil Function modification
- Phase and Magnitude Aberrations via Zernike expansion
- Any PSF can be converted to an interpolated PSF for faster generation at new positions   

## Design

The high-level interface is designed to facilitate generation of synthetic data as would be seen by an Array Detector (e.g. Camera).  The PSF is considered a probability distribution that is normalized across 2D sections.  Calculating the PSF at a location follows the convention from  [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) where a distribution is created, and the PDF is calculated at a location.  

## Examples

### Airy PSF 

```julia
using MicroscopePSFs
PSF=MicroscopePSFs

# Create an Airy PSF
na=1.2
λ=.6 
pixelsize=.1
p=PSF.Airy2D(na,λ,pixelsize)

# calculate the PSF at a point
PSF.pdf(p,(0,0),(0,0))

# calculate the PSF in a region
sz=16
roi=[(i,j) for i=1:sz, j=1:sz]
im1=PSF.pdf(p,roi,(sz/2,sz/2))
```

### 3D Scalar PSF with Astigmatims

```julia
using MicroscopePSFs
PSF=MicroscopePSFs

# Zernike Magnitude and Phase Coefficients 
z_mag=[1.0]
z_phase=zeros(10)
z_phase[6]=1 # astigmatism
z=PSF.ZernikeCoefficients(z_mag,z_phase)

# Create a scalar PSF
na=1.2
n=1.3
λ=.6 
pixelsize=.1

p=PSF.Scalar3D(na,λ,n,pixelsize;z=z)

# calculate the PSF in a region
sz=32
roi=[(i,j,0) for i=-sz/2:(sz/2-1), 
    j=-sz/2:(sz/2-1)] 
im=PSF.pdfₐ(p,roi,(0.0,0.0,-10.0))
```


## Future Development

*MicroscopePSFs* should be considered under development and the interface may change as we build the JuliaSMLM ecosystem.  Comments and Feature requests are welcome via an issue.  

### Future Features
- Vector PSFs
- Super-critical angle PSFs
- Integration across finite pixels sizes using sub-sampling.  
