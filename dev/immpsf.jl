using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots
using FisherInfoPoisson

#filename = raw"Y:\Projects\Super Critical Angle Localization Microscopy\ZStack_TIRF_01-10-23\ZStack_psfmodel_zernike_vector_single.h5"


# define PSF parameters
n = [1.33,1.52,1.52]
na = 1.49
λ = 0.69
pixelsize = 0.05
zstage = 0.0 # stage position um
mag=[1.0]
phase=zeros(10)
phase[6]=0.0 # astigmatism
z=PSF.ZernikeCoefficients(mag,phase)

p = PSF.ImmPSF(na, λ, n, pixelsize; zstage=zstage, ksize=256,mvtype="stage")
#p,_,_,_ = PSF.importpsf(filename,"immPSF",zstage=zstage,mvtype="stage")  

h1 = p.pupilfunction[1]
p1 = heatmap(h1.pupil[:,:,1], aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:grays)
p2 = heatmap(h1.pupil[:,:,2], aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:grays)
plot(p1,p2,layout=(1,2))


sz = 20
roi=[(x,y,0) for x=-sz/2:(sz/2-1), y=-sz/2:(sz/2-1)] 
pos = (0.0,0.0,0.0)
im1=PSF.pdf(p,roi,pos)
p1 = heatmap(im1, aspectratio=:equal, yflip = true, colorbar=:none,axis = nothing,showaxis=false)
h = plot(p1,size=(400,400))
savefig(h,"freedipole_psf.png")

# Generate a PSF stack
sz = 20
roi=[(x,y,k) for x=0:sz-1,y=0:sz-1,k=0:0]
xe = sz/2
ye = sz/2
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=0.5:-0.05:-0.5]


for j=eachindex(pos)
    im1=PSF.pdf(p,roi,pos[j])
    plt=heatmap(im1[:,:,1], aspectratio=:equal, yflip = true)
    zpos = pos[j][3]
    plot!(plt,title="PSF, z: $zpos")
    display(plt)
    sleep(.01)
    print(sum(im1))
    print("\n")
end

# calculate CRLB    
model(θ, (x, y, z)) = θ[4] * PSF.pdf(p, (x, y, z), (θ[1], θ[2], θ[3])) + θ[5]
θ = [pos[1]...,1000,10]    
fi,crlb=FisherInfoPoisson.calcFI(model,θ,roi,prior=Inf)



