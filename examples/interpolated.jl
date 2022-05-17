# Make an interpolated PSF 

using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots

# Microscope
na=1.2
n=1.3
λ=.6 
pixelsize=.1
sz=16 
roi=[(i,j) for i=1:sz,j=1:sz]

##  Airy2D PSF
p=PSF.Airy2D(na,λ,pixelsize)
ip=PSF.InterpolatedPSF(p,(sz*2,sz*2);subsampling=2)
im=PSF.pdf(ip,roi,(sz/2,sz/2))
heatmap((im))
sum(im) #check normalization 

@time im=PSF.pdf(p,roi,(sz/2,sz/2));
@time im=PSF.pdf(ip,roi,(sz/2,sz/2));
@time PSF.pdf!(im,ip,roi,(sz/2,sz/2));

## Scalar3D PSF 
mag=[1.0]
phase=zeros(10)
phase[6]=1 # astigmatism
z=PSF.ZernikeCoefficients(mag,phase)
roi=[(i,j,k) for i=1:sz,j=1:sz,k=0:0]
p=PSF.Scalar3D(na,λ,n,pixelsize;z)
ip=PSF.InterpolatedPSF(p,(sz*2,sz*2,1);subsampling=2)
@time im=PSF.pdf(ip,roi,(sz/2,sz/2,0))
plt=heatmap((im[:,:,1]))
display(plt)

#check normalization 
display(sum(im))

# Look at psfₐ
psfₐ=PSF.pdfₐ(ip,roi,(sz/2,sz/2,0))
plt2=heatmap(abs2.(psfₐ[:,:,1])) 
display(plt2)

#compare to direct psf calc 
plt3=heatmap(abs2.(psfₐ[:,:,1]).-im[:,:,1]) 


























