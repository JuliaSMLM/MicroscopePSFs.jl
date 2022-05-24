# A few tests and examples of zernike polynomials

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots 

rsz=128
ρ=[ sqrt(ii^2+jj^2)/rsz for ii=-rsz:rsz-1, jj=-rsz:rsz-1]
ϕ=[ atan(ii,jj) for ii=-rsz:rsz-1, jj=-rsz:rsz-1]

## 
n=5
l=-3

PSF.nl2osa(n,l)

im=PSF.zernikepolynomial.(n,l,ρ,ϕ)
heatmap(im)

## OSA index
for j=0:20
    im=PSF.zernikepolynomial.(j,ρ,ϕ)
    plt=heatmap(im)
    n,l=PSF.osa2nl(j)
    plot!(plt,title="OSA Index: $j, n: $n l:$l")
    display(plt)
    sleep(.1)
end


