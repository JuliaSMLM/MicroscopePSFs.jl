# Create and show a PSF with Zernike magnitude and phase abbeerations 

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

mag=[1.0]
phase=zeros(10)
phase[8]=0.5 # astigmatism

z=PSF.ZernikeCoefficients(mag,phase)

# Create a scalar PSF
na=1.35
n=1.406
λ=.69 
pixelsize=.13

p=PSF.Scalar3D(na,λ,n,pixelsize;z=z)

# Or Create a 2D Gaussian approximation
#p=PSF.Gauss2D(p)

# calculate the PSF at a point
PSF.pdf(p,(0,0,0),(0.0,0.0,0.0))

# calculate the PSF in a region
sz=32
roi=[(x,y,0) for y=-sz/2:(sz/2-1), 
    x=-sz/2:(sz/2-1)] 

# out of focus by up to 0.5 physical unit (usually micron)
zrange=cat(dims=1,collect(LinRange(-0.5,0.5,5))
            ,collect(LinRange(0.5,-0.5,5)))
anim = @animate for z ∈ zrange
   heatmap(PSF.pdf(p,roi,(0.0,0.0,z)), aspectratio=:equal, 
   yflip = true, colorbar=:none)
end
gif(anim, fps = 1)


