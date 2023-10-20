using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots


# Create a scalar PSF
na=1.49
n=[1.33,1.52,1.52] # refractive indices (sample medium, cover glass, immersion)
λ=.69
pixelsize=.106
dipole_ang = [0,0].*pi./180
p=PSF.Dipole3D(na,λ,n,pixelsize,dipole_ang;ksize=256)

# look at pupil
h1 = p.pupilfunctionx.pupil[:,:,1].^2+p.pupilfunctiony.pupil[:,:,1].^2
p1 = heatmap(h1, aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:viridis)

h1 = (p.pupilfunctionx.pupil[:,:,1].*p.apodization[:,:,1]).^2+(p.pupilfunctiony.pupil[:,:,1].*p.apodization[:,:,1]).^2
h1 = p.apodization[:,:,1].^2
p1 = heatmap(h1, aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:viridis)


h1 = p.pupilfunctionx
p1 = heatmap(h1.pupil[:,:,1], aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:grays)
p2 = heatmap(h1.pupil[:,:,2], aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:grays)
plot(p1,p2,layout=(1,2))

#calculate the PSF in a region
sz=16
roi=[(x,y,0) for x=-sz/2:(sz/2-1), y=-sz/2:(sz/2-1)] 

pos_emitter = (0.0,0.0,0.0)
p.electricfield = 'x'
imx=PSF.pdf(p,roi,pos_emitter)
p.electricfield = 'y'
imy=PSF.pdf(p,roi,pos_emitter)

#look at psf
heatmap((imx+imy), aspectratio=:equal, yflip = true, colorbar=:none)


xe = 0.0
ye = 0.0
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=-0.5:0.1:0.5]
out = zeros(sz,sz,length(pos))
for j=eachindex(pos)
    p.electricfield = 'x'
    im=PSF.pdf(p,roi,pos[j])
    plt=heatmap(im[:,:,1], aspectratio=:equal, yflip = true)
    zpos = pos[j][3]
    plot!(plt,title="PSF, z: $zpos")
    display(plt)
    sleep(.1)
    print(sum(im))
    print("\n")
end