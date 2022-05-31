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

#calculate the PSF in a region
sz=16
roi=[(x,y,0) for x=-sz/2:(sz/2-1), 
    y=-sz/2:(sz/2-1)] 

im=PSF.pdf(p,roi,(0.0,0.0,0.0))

#look at psf
plt=heatmap((im), aspectratio=:equal, yflip = true)
display(plt)

#check normalization 
sum(im)


# out of focus by up to 1 physical unit (usually micron)
zrange=collect(LinRange(-1,1,5))
for z ∈ zrange
   display(heatmap(PSF.pdf(p,roi,(0.0,0.0,z))))
end





