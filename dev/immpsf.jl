using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots
using FisherInfoPoisson

filename = raw"Y:\Personal Folders\Sheng\data\PSF Engineering Microscope\PR_10-04-2022_UAF_7px_Shift_200nm_beads_Gain0_PupilSeq\ZStack_-2to+2um_PupilSize125psfmodel_LL_pupil_vector_single.h5"

zstage = 0.0 # stage position 1um

p = PSF.importpsf(filename,"immPSF",zstage=zstage)  


h1 = p.pupilfunction[5]
p1 = heatmap(h1.pupil[:,:,1], aspectratio=:equal, yflip = true, axis = nothing,showaxis=false)
p2 = heatmap(h1.pupil[:,:,2], aspectratio=:equal, yflip = true, axis = nothing,showaxis=false)
plot(p1,p2,layout=(1,2))


# Generate a PSF stack
sz = 21
roi=[(x,y,k) for x=0:sz-1,y=0:sz-1,k=0:0]
xe = sz/2
ye = sz/2
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=0:0.2:1]


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



