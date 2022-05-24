# Create and show an amplitude PSF with Zernike magnitude and phase abbeerations 

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

z_mag=[1.0]
z_phase=zeros(10)
z_phase[6]=1 # astigmatism

z=PSF.ZernikeCoefficients(z_mag,z_phase)

# Create a scalar PSF
na=1.2
n=1.3
λ=.6 
pixelsize=.1

p=PSF.Scalar3D(na,λ,n,pixelsize;z=z)

# Or Create a 2D Gaussian approximation
#p=PSF.Gauss2D(p)

# calculate the PSF at a point
PSF.pdfₐ(p,(0,0,0),(0.0,0.0,0.0))

# calculate the PSF in a region
sz=32
roi=[(i,j,0) for i=-sz/2:(sz/2-1), 
    j=-sz/2:(sz/2-1)] 

# check normalization    
im=PSF.pdfₐ(p,roi,(0.0,0.0,0.0))
println("integrated intensity: ", sum(abs2.(im)))

# look at psf. Note x,y are pixels, z is physical unit. 
im=PSF.pdfₐ(p,roi,(0.0,0.0,-1.0))
display(heatmap(abs.(im)))
display(heatmap(angle.(im)))


