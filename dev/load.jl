
using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using CairoMakie
CM = CairoMakie


filename = raw"/mnt/nas/lidkelab/Projects/TIRF Demo/fernando/06-03-24/psfmodel_zernike_vector_single_smallz.h5"
p, PSFstack, z, h = PSF.importpsf(filename,"splinePSF")  
#p,_,_,_ = PSF.importpsf(filename,"immPSF",zstage=1.0,mvtype="stage")  

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
roi=[(y,x,k) for y=0:sz-1,x=0:sz-1,k=0:0]
xe = sz/2
ye = sz/2
pos = [(y,x,k) for y=ye:ye,x=xe:xe,k=-0.4:0.1:0.4]

#@enter im=PSF.pdfₐ(ip,roi[1],pos[1])

for j=eachindex(pos)
    ims=PSF.pdf(ip,roi,pos[j])
    zpos = pos[j][3]
    fig = Figure(size = (300, 300))
    ax = CM.Axis(fig[1, 1],title="PSF, z: $zpos",aspect=1,xreversed=true)
    hm = CM.heatmap!(ax, ims[:,:,1],colormap=:inferno)
    #hidedecorations!(ax)
    fig
    display(fig)
    sleep(.1)
    print(sum(ims))
    print("\n")
end



# load FD PSF

filename = raw"Y:\Projects\Super Critical Angle Localization Microscopy\Data\10-06-2023\Data4\psf_kmed_tilt2_insitu_zernike_single_5D_tilt.h5"
p, PSFstack, z, h = PSF.importpsf(filename,"splinePSF_FD");  

psffile = splitext(filename)[1]*"_1.jld2"
PSF.save(psffile,p)
ip = PSF.load(psffile);

sz = 20
roi=[(x,y,k) for x=0:sz-1,y=0:sz-1,k=0:0]
xe = sz/2
ye = sz/2
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=0:0.1:1.3]
cor = (245,0)


for j=eachindex(pos)
    ims=PSF.pdf(ip,roi,pos[j],cor)
    zpos = pos[j][3]
    fig = Figure(size = (100, 100))
    ax = CM.Axis(fig[1, 1],title="PSF, z: $zpos",aspect=1,xreversed=true)
    hm = CM.heatmap!(ax, ims[:,:,1],colormap=:inferno)
    #hidedecorations!(ax)
    fig
    display(fig)
    sleep(.1)
    print(sum(ims))
    print("\n")
end

# load stage tilt 4D PSF
filename = raw"Y:\Projects\Super Critical Angle Localization Microscopy\Data\10-06-2023\Data4\psf_kmed2_insitu_zernike_single_4D_tilt.h5"
p, PSFstack, z, h = PSF.importpsf(filename,"splinePSF_tilt");  

psffile = splitext(filename)[1]*"_1.jld2"
PSF.save(psffile,p)
ip = PSF.load(psffile);

sz = 20
roi=[(y,x,k) for y=0:sz-1,x=0:sz-1,k=0:0]
xe = sz/2
ye = sz/2
pos = [(y,x,k) for y=ye:ye,x=xe:xe,k=0:0.02:0.2]
stpos = 0.0




for j=eachindex(pos)
    ims=PSF.pdf(ip,roi,pos[j],stpos)
    zpos = pos[j][3]
    fig = CM.Figure(size = (100, 100))
    ax = CM.Axis(fig[1, 1],title="PSF, z: $zpos",aspect=1,yreversed=true)
    hm = CM.heatmap!(ax, ims[:,:,1]',colormap=:inferno)
    #hidedecorations!(ax)
    fig
    display(fig)
    sleep(.1)
    print(sum(ims))
    print("\n")
end


img = PSFstack[:,:,100,1]
fig = Figure(size = (100, 100))
ax = CM.Axis(fig[1, 1],aspect=1,yreversed=true)
hm = CM.heatmap!(ax, img',colormap=:inferno)
fig