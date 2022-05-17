using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots


# Create a scalar PSF
na=1.2
n=1.3
λ=.6 
pixelsize=.1

p=PSF.Scalar3D(na,λ,n,pixelsize)

# Or Create a 2D Gaussian approximation
#p=PSF.Gauss2D(p)

#calculate the PSF at a point
PSF.pdf(p,(0,0,0),(0.0,0.0,0.0))

#calculate the PSF in a region
sz=16
roi=[(i,j,0) for i=-sz/2:(sz/2-1), 
    j=-sz/2:(sz/2-1)] 

im=PSF.pdf(p,roi,(0.0,0.0,0.0))

#look at psf
plt=heatmap((im))
display(plt)
#check normalization 
sum(im)


# out of focus
im=PSF.pdf(p,roi,(0.0,0.0,10.0))

#look at psf
heatmap((im))




