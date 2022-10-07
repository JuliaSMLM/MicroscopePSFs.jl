
using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots


filename = raw"examples\psfmodel_tetropod_zernike_single.h5"
p,PSFstack = PSF.importpsf(filename)  
 

# interpolate psf: linear
#sz = 16
#ip=PSF.InterpolatedPSF(p,(sz*2,sz*2,1);subsampling=1) # z unit is in micron

# interpolate psf: cubic
ip = PSF.SplinePSF(PSFstack) # z unit is in pixel

# Generate a PSF stack
sz = 16
roi=[(x,y,k) for x=0:sz-1,y=0:sz-1,k=0:0]
xe = 8
ye = 8
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=-1:0.5:1]


for j=eachindex(pos)
    im=PSF.pdf(ip,roi,pos[j])
    plt=heatmap(im[:,:,1], aspectratio=:equal, yflip = true)
    zpos = pos[j][3]
    plot!(plt,title="PSF, z: $zpos")
    display(plt)
    sleep(.1)
    print(sum(im))
    print("\n")
end







