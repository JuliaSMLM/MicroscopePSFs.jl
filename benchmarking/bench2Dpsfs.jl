# Make an interpolated PSF 

using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots
using BenchmarkTools
using ProfileView

# Microscope
na=1.2
n=1.3
λ=.6 
pixelsize=.1
sz=16
roi=[(i,j) for i=1:sz,j=1:sz]

##
function speedtest(p,roi,sz)
    for ii=1:1000
        image=PSF.pdf(p,roi,(sz/2,sz/2))
    end
    return image
end

##  Airy2D PSF
p=PSF.Airy2D(na,λ,pixelsize)
@benchmark image=PSF.pdf(p,roi,(sz/2,sz/2))
image=PSF.pdf(p,roi,(sz/2,sz/2));
@benchmark PSF.pdf!(image,p,roi,(sz/2,sz/2))
@btime speedtest(p,roi,sz)
ProfileView.@profview  speedtest(p,roi,sz)

## Gauss2D
p=PSF.Gauss2D(p)
@benchmark image=PSF.pdf(p,roi,(sz/2,sz/2))
image=PSF.pdf(p,roi,(sz/2,sz/2));
@benchmark PSF.pdf!(image,p,roi,(sz/2,sz/2))
@btime speedtest(p,roi,sz)
ProfileView.@profview  speedtest(p,roi,sz)

## Interpolated 
ip=PSF.InterpolatedPSF(p,(sz*2,sz*2))
@benchmark image=PSF.pdf(ip,roi,(sz/2,sz/2))
image=PSF.pdf(ip,roi,(sz/2,sz/2))
@benchmark PSF.pdf!(image,ip,roi,(sz/2,sz/2))
ProfileView.@profview  speedtest(ip,roi,sz)

