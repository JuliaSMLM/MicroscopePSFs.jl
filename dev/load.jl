
using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots


filename = raw"Y:\Projects\Super Critical Angle Localization Microscopy\ZStack_TIRF_01-10-23\ZStack_psfmodel_zernike_vector_single.h5"
p,PSFstack,z = PSF.importpsf(filename,"splinePSF",zstage = 1.0)  
 

psffile = splitext(filename)[1]*".jld2"

PSF.save(psffile,p)
ip = PSF.load(psffile)

# interpolate psf: linear
#sz = 16
#ip=PSF.InterpolatedPSF(p,(sz*2,sz*2,1);subsampling=1) # z unit is in micron

# interpolate psf: cubic
#ip = PSF.SplinePSF(PSFstack) # z unit is in pixel

# Generate a PSF stack
sz = 20
roi=[(x,y,k) for x=0:sz-1,y=0:sz-1,k=0:0]
xe = sz/2
ye = sz/2
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=0:0.1:2]

#@enter im=PSF.pdfₐ(ip,roi[1],pos[1])
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






