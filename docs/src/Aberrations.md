# Aberrations

*MicroscopePSFs* includes a class of PSFs that are generated via integration over a complex pupil function, $P(k_x,k_y)$. Aberrations can be included in this calculation by modification of the pupil function.  This can be done by direct modification of the pupil or by expansion into Zernike modes.  

## The Pupil Function 




## Zernike Expansion

Pupil function based PSFs can be created with a pupil magnitude and phase that are each given by a sum of Zernike polynomials.  These are specified by a [`MicroscopePSFs.ZernikeCoefficients`](@ref) structure that is passed to the consructor:


The fields of this structure hold the coefficients of Zernike expansion using the [OSA/ANSI](https://en.wikipedia.org/wiki/Zernike_polynomials#OSA/ANSI_standard_indices) linear index.  These vectors are one based where e.g. `mag[1]` holds the coefficient for the `j = 0` linear index.  

Zernike polynomials are normalized such that coefficients are equal to the root-mean-square error across the unit disk.  This normalization is:

$\int_0^{2\pi}\int_0^1 Z^2  \rho d\rho d\phi = \pi$


### Example 

A Tetrapod type PSF using a mixture of 1st and 2nd order astigmatism. 
```@setup
using Pkg
Pkg.add("Plots")
```

```@example
using Plots
using MicroscopePSFs

mag=[1.0]
phase=zeros(14)
phase[6]=.5 #osa index 5
phase[14]=-.5 #osa index 13
z=MicroscopePSFs.ZernikeCoefficients(mag,phase)

na=1.2
n=1.3
λ=.6 
pixelsize=.1

p=MicroscopePSFs.Scalar3D(na,λ,n,pixelsize;z=z)

sz=32  # hide
roi=[(x,y,0) for y=-sz/2:(sz/2-1), x=-sz/2:(sz/2-1)]  # hide

zrange=cat(dims=1,collect(LinRange(-1,1,9)),collect(LinRange(1,-1,9)))  # hide
anim = @animate for z ∈ zrange  # hide
   heatmap(MicroscopePSFs.pdf(p,roi,(0.0,0.0,z)), aspectratio=:equal, yflip = true, colorbar=:none)  # hide
end  # hide
gif(anim, "tetrapod.gif",fps = 5) # hide
```


