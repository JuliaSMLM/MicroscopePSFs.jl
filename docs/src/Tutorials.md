# Import, save and load PSFs

## Import

The `importpsf` function imports the PSF structure from the [uiPSF](https://github.com/ries-lab/uiPSF) package as a `MicroscopePSFs` type. This function requires inputs of the PSF filename and the PSF type. 

```@docs
MicroscopePSFs.importpsf
```
### Import Example:
```julia
using MicroscopePSFs
PSF=MicroscopePSFs

filename = "psf_from_uiPSF.h5"
p,PSFstack,z,h = PSF.importpsf(filename,"splinePSF") 
p,PSFstack,z,h = PSF.importpsf(filename,"immPSF",zstage=1.0,mvtype="bead")
p,PSFstack,z,h = PSF.importpsf(filename,"scalar3D") 
```
## Save
Save the `MicroscopePSFs` type as a .jld2 file. The example below first creates a `Scalar3D` PSF type and then save the PSF as a .jld2 file. 

```@docs
MicroscopePSFs.save
```
### Save Example
```julia
# Create a scalar PSF
na=1.2          # Numerical Aperture
n=1.3           # Refractive Index
λ=.6            # Wavelength (micron)
pixelsize=.1    # Pixel Size (micron)

p=PSF.Scalar3D(na,λ,n,pixelsize)

# Save it
PSF.save("psf.jld2",p)
```
## Load

Load the `MicroscopePSFs` type from a .jld2 file.

```@docs
MicroscopePSFs.load
```
### Load Example
```julia
p = PSF.load("psf.jld2")
```