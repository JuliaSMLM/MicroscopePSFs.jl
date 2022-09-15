
using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots


filename = raw"examples\psfmodel_test_insitu_single.h5"
p,PSFstack = PSF.loadh5(filename)  


# interpolate psf: linear
#sz = 2
#ip=PSF.InterpolatedPSF(p,(sz*2,sz*2,1);subsampling=2)

# interpolate psf: cubic
ip = PSF.SplinePSF(PSFstack)

# Generate a PSF stack
sz = 11
roi=[(x,y,k) for x=1:sz,y=1:sz,k=0:0]
xe = 6.2
ye = 5.5
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=-20:0.5:5]

out = zeros(sz,sz,length(pos))

for ii=eachindex(pos)
    out[:,:,ii]=PSF.pdf(ip,roi,pos[ii])
end
heatmap(out, aspectratio=:equal, yflip = true, colorbar=:none)



heatmap(out[:,:,2], aspectratio=:equal, yflip = true, colorbar=:none)







