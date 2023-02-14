## Create a psf, calculate in a ROI and display
using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

# Create an Airy PSF
na=1.2
λ=.6 
pixelsize=.1
p_airy=PSF.Airy2D(na,λ,pixelsize)

# Or Create a 2D Gaussian approximation
p_gauss=PSF.Gauss2D(p_airy)

#pick one
p=p_gauss
#p=p_airy

# calculate the PSF at a point
PSF.pdf(p,(0,0),(0,0))

# calculate the PSF in a region
sz=16
roi=[(x,y) for y=1:sz, x=1:sz]

im1=PSF.pdf(p,roi,(sz/2+3,sz/2))
# check normalization 
sum(im1)

# look at psf
plt=heatmap(im1, aspectratio=:equal, yflip = true)
display(plt)

# calculate for multiple emitters
x_emitters=[(sz/2-2,sz/2),(sz/2+2,sz/2)]
im2=PSF.pdf(p,roi,x_emitters)
heatmap(im2,aspectratio=:equal, yflip = true)








