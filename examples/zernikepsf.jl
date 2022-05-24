# Create and show a PSF with Zernike magnitude and phase abbeerations 

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

mag=[1.0]
phase=zeros(10)
phase[6]=1 # astigmatism

z=PSF.ZernikeCoefficients(mag,phase)

# Create a scalar PSF
na=1.2
n=1.3
λ=.6 
pixelsize=.1

p=PSF.Scalar3D(na,λ,n,pixelsize;z=z)

# Or Create a 2D Gaussian approximation
#p=PSF.Gauss2D(p)

# calculate the PSF at a point
PSF.pdf(p,(0,0,0),(0.0,0.0,0.0))

# calculate the PSF in a region
sz=32
roi=[(i,j,0) for i=-sz/2:(sz/2-1), 
    j=-sz/2:(sz/2-1)] 

# out of focus by up to 1 physical unit (usually micron)
zrange=collect(LinRange(-1,1,9))
for z ∈ zrange
   display(heatmap(PSF.pdf(p,roi,(0.0,0.0,z))))
end


