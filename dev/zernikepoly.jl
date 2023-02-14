# A few tests and examples of zernike polynomials

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots 
using Statistics

rsz=128
ρ=[ sqrt(ii^2+jj^2)/rsz for ii=-rsz:rsz-1, jj=-rsz:rsz-1]
ϕ=[ atan(ii,jj) for ii=-rsz:rsz-1, jj=-rsz:rsz-1]

## 
n=3
l=-1

PSF.nl2noll(n,l)

im=PSF.zernikepolynomial.(n,l,ρ,ϕ)
heatmap(im,aspectratio=:equal, yflip = true)


## OSA index
for j=1:21
    im=PSF.zernikepolynomial.(j,ρ,ϕ,linearindex="Noll")
    plt=heatmap(im, aspectratio=:equal, yflip = true)
    n,l=PSF.noll2nl(j)
    plot!(plt,title="OSA Index: $j, n: $n l:$l")
    display(plt)
    sleep(.1)
    print("\n")
    display(sum(abs2.(im)) / (rsz)^2 )
end


