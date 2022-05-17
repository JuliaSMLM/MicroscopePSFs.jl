
using Revise
using MicroscopePSFs
using BenchmarkTools
using ProfileView
using LoopVectorization
PSF=MicroscopePSFs

# Create a scalar PSF
na=1.2
n=1.3
λ=.6 
pixelsize=.1

p=PSF.Scalar3D(na,λ,n,pixelsize)

sz=8
roi=[(i,j,k) for i=-sz/2*pixelsize:pixelsize:pixelsize*(sz/2-1), 
    j=-sz/2*pixelsize:pixelsize:pixelsize*(sz/2-1),k=-1:.1:1]

@btime im=PSF.pdf.(p,roi)

ProfileView.@profview PSF.pdf.(p,roi)
ProfileView.@profview PSF.pdf.(p,roi)

# interpolated 







