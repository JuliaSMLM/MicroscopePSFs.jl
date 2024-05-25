using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using Plots
using CairoMakie
CM = CairoMakie
# Create a scalar PSF
na=1.49
n=[1.33,1.52,1.52] # refractive indices (sample medium, cover glass, immersion)
λ=.69
pixelsize=.1
dipole_ang = [0,0].*pi./180
p=PSF.Dipole3D(na,λ,n,pixelsize,dipole_ang;ksize=128,excitationfield=1.0)

# look at pupil
#h1 = p.pupilfunctionx.pupil[:,:,1].^2+p.pupilfunctiony.pupil[:,:,1].^2
#p1 = heatmap(h1, aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:viridis)
#
#h1 = (p.pupilfunctionx.pupil[:,:,1].*p.apodization[:,:,1]).^2+(p.pupilfunctiony.pupil[:,:,1].*p.apodization[:,:,1]).^2
#h1 = p.apodization[:,:,1].^2
#p1 = heatmap(h1, aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,c=:viridis)
hx = p.pupilfunctionx.pupil[:,:,1].*exp.(im*p.pupilfunctionx.pupil[:,:,2])
hy = p.pupilfunctiony.pupil[:,:,1].*exp.(im*p.pupilfunctiony.pupil[:,:,2])
normf = sqrt(real(sum(hx.*conj(hx))+sum(hy.*conj(hy))))



dipole_ang = [90,0].*pi./180
p=PSF.Dipole3D(na,λ,n,pixelsize,dipole_ang;normf =normf, ksize=128)
hx = p.pupilfunctionx.pupil[:,:,1].*exp.(im*p.pupilfunctionx.pupil[:,:,2])
hy = p.pupilfunctiony.pupil[:,:,1].*exp.(im*p.pupilfunctiony.pupil[:,:,2])
fig = Figure(size = (600, 300))
ax = CM.Axis(fig[1, 1],title="Ex",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, abs.(hx),colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 2],title="Ey",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, abs.(hy),colormap=:inferno)
hidedecorations!(ax)
fig


#p1 = heatmap(abs.(p.pupilfunctionx.pupil[:,:,1]), aspectratio=:equal, yflip = true, axis = nothing,showaxis=false, colorbar=false,c=:inferno)
#p2 = heatmap(abs.(p.pupilfunctiony.pupil[:,:,1]), aspectratio=:equal, yflip = true, axis = nothing,showaxis=false,colorbar=false,c=:inferno)
#h = plot(p1,p2,layout=(1,2),size=(800,400))
#savefig(h,"dipole_pupil_90.png")

#calculate the PSF in a region
sz=40 # there is a permutation in PSF calculation
roi=[(x,y,0) for x=-sz/2:(sz/2-1), y=-sz/2:(sz/2-1)] 

pos_emitter = (0.0,0.0,0.0)
p.electricfield = 'x'
imx=PSF.pdf(p,roi,pos_emitter)
p.electricfield = 'y'
imy=PSF.pdf(p,roi,pos_emitter)

#look at psf
fig = Figure(size = (800, 400))
ax = CM.Axis(fig[1, 1],title="Ex psf",aspect=1)
hm = CM.heatmap!(ax, imx,colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 2],title="Ey psf",aspect=1)
hm = CM.heatmap!(ax, imy,colormap=:inferno)
hidedecorations!(ax)
fig


xe = 0.0
ye = 0.0
pos = [(x,y,k) for x=xe:xe,y=ye:ye,k=0.0:0.1:0.5]
out = zeros(sz,sz,length(pos))
for j=eachindex(pos)
    p.electricfield = 'x'
    imx=PSF.pdf(p,roi,pos[j])
    p.electricfield = 'y'
    imy=PSF.pdf(p,roi,pos[j])
    ims = imx.+imy
    zpos = pos[j][3]

    fig = Figure(size = (400, 400))
    ax = CM.Axis(fig[1, 1],title="PSF, z: $zpos",aspect=1)
    hm = CM.heatmap!(ax, ims,colormap=:inferno)
    hidedecorations!(ax)
    fig
    display(fig)
    sleep(.1)
    print(sum(ims))
    print("\n")
end


ang = []
Nθ = 12
θ = Array(range(0, pi, Nθ))
dθ = θ[2] - θ[1]
for j in eachindex(θ)
    dϕ = dθ / sin(θ[j])
    Nϕ = round(Int, 2 * pi / dϕ)
    for i in 0:Nϕ
        if i * dϕ < (2 * pi - dϕ)
            append!(ang, [(θ[j], i * dϕ)])
        end
    end
end
append!(ang, [(0.0, 0.0)])
append!(ang, [(pi, 0.0)])

ang_array = [tup[j] for tup in ang,j = 1:2]
α = ang_array[:,1]
β = ang_array[:,2]
z = cos.(α)    
x = sin.(α).*cos.(β)
y = sin.(α).*sin.(β)
r = hcat(x,y,z)
r1 = sortslices(r,dims=1,by=x->x[3],rev=true)
fig = Figure(size = (400, 400))
ax = CM.Axis3(fig[1, 1],aspect = :equal,xlabel="x (nm)",ylabel="y (nm)",zlabel="z (nm)")
CM.scatter!(ax,r1[:,1],r1[:,2],r1[:,3],markersize=10,marker=:circle)
CM.lines!(ax,r1[:,1],r1[:,2],r1[:,3])
fig

