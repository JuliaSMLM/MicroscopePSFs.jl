## Create a psf, calculate in a ROI and display

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

# Create an Airy PSF
na=1.2
λ=.6 
pixelsize=.1
p=PSF.Airy2D(na,λ,pixelsize)

# Or Create a 2D Gaussian approximation
p=PSF.Gauss2D(p)

# calculate the PSF at a point
PSF.pdf(p,(0,0),(0,0))

# calculate the PSF in a region
sz=16
roi=[(i,j) for i=1:sz, j=1:sz]

im1=PSF.pdf(p,roi,(sz/2,sz/2))
# check normalization 
sum(im1)

# look at psf
plt=heatmap(im1)
display(plt)

# calculate for multiple emitters
x_emitters=[(sz/2-2,sz/2),(sz/2+2,sz/2)]
im2=PSF.pdf(p,roi,x_emitters)
heatmap(im2)








