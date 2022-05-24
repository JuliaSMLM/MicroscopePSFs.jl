# Make an interpolated PSF 

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

# Microscope
na=1.2
n=1.3
λ=.6 
pixelsize=.1
sz=16 

## Scalar3D PSF 
mag=[1.0]
phase=zeros(10)
phase[6]=1 # astigmatism
z=PSF.ZernikeCoefficients(mag,phase)
p=PSF.Scalar3D(na,λ,n,pixelsize;z)
ip=PSF.InterpolatedPSF(p,(sz*2,sz*2,.2);subsampling=2)

roi=[(i,j,k) for i=1:sz,j=1:sz,k=0:0]

@time im=PSF.pdf(p,roi,(sz/2,sz/2,.2))
@time imp=PSF.pdf(ip,roi,(sz/2,sz/2,.2))

# check normalization 
display(sum(im))

#  display both
plt=heatmap(cat(dims=2,im[:,:,1],imp[:,:,1]); aspectratio = :equal, yflip = true)
display(plt)




























