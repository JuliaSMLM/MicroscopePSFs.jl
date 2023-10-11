# Interpolated PSFs

## Overview 
Interpolated PSFs can be usefull when the same PSF is needed to generate models with various different source locations. 
Interpolated PSFs use the [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) package. While the Airy and Gauss PSFs
can be interpolated, speed benifits primarily come from avoiding the integral over the pupil function in the pupil function based PSFs 
such as `Scaler3D`.  Due to sub-sampling, 3D PSFs can be slow to generate.    


[`MicroscopePSFs.InterpolatedPSF`](@ref)


## Examples

```julia
using MicroscopePSFs
PSF=MicroscopePSFs

na=1.2
n=1.3
λ=.6 
pixelsize=.1
sz=16 

# Create the PSF 
p=PSF.Scalar3D(na,λ,n,pixelsize)

#Build Interpolation
maxrange=(sz*2,sz*2,1)
ip=PSF.InterpolatedPSF(p,maxrange)

im=PSF.pdf(ip,roi,(sz/2,sz/2,0))
```


