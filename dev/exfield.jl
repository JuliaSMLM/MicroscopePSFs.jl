using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using CairoMakie
CM = CairoMakie
using ImageView
using GLMakie
GM = GLMakie

na=1.35
n=[1.406,1.52,1.33] # refractive indices (immersion , cover glass, sample medium)
λ=.642
pixelsize=.04

# Focusfield(nₐ, λ, n::Vector, pixelsize;
#     normf=1.0, 
#     ksize=256,
#     electricfield='x',
#     excitationfield=[1.0,0],
#     fθ = (x)->0.0, 
#     fϕ = (x)->0.0)

p=PSF.Focusfield(na,λ,n,pixelsize; ksize=128,excitationfield=[1.0,1*im]) 

function phase_theta(x, na, nimm,ringsz)
    θmax = asin(na/nimm)
    Nring = length(ringsz)
    if x < ringsz[1]*θmax
        return π
    end

    for j in 2:Nring
        if (x > ringsz[j-1]*θmax) && (x < ringsz[j]*θmax)
            return (j-2)*π
        end
    end
end

function phase_phi(x)
    return x
end

fθ = (x) -> phase_theta(x, na, n[1], [0.4, 1.0])

p=PSF.Focusfield(na,λ,n,pixelsize; ksize=128,excitationfield=[1.0,1*im],fθ=fθ) 

function focus_gui(na,λ,n,pixelsize)
    imx_xy = Observable{Any}(0.0)
    imy_xy = Observable{Any}(0.0)
    imz_xy = Observable{Any}(0.0)
    imall_xy = Observable{Any}(0.0)
    I_x = Observable{Any}(0.0)

    imx_xz = Observable{Any}(0.0)
    imy_xz = Observable{Any}(0.0)
    imz_xz = Observable{Any}(0.0)
    imall_xz = Observable{Any}(0.0)
    I_z = Observable{Any}(0.0)

    fig = GM.Figure(resolution = (700,400), title = "focus field")
    ax_Ex_xz = GM.Axis(fig,aspect = 1,title = "Ex",yreversed=true)
    ax_Ey_xz = GM.Axis(fig,aspect = 1,title = "Ey",yreversed=true)
    ax_Ez_xz = GM.Axis(fig,aspect = 1,title = "Ez",yreversed=true)
    ax_Etot_xz = GM.Axis(fig,aspect = 1,title = "total",yreversed=true)
    ax_profile_z = GM.Axis(fig,title = "Iz")
    fig[1,:] = hgrid!(ax_Ex_xz,ax_Ey_xz,ax_Ez_xz,ax_Etot_xz,ax_profile_z)

    ax_Ex_xy = GM.Axis(fig,aspect = 1,yreversed=true)
    ax_Ey_xy = GM.Axis(fig,aspect = 1,yreversed=true)
    ax_Ez_xy = GM.Axis(fig,aspect = 1,yreversed=true)
    ax_Etot_xy = GM.Axis(fig,aspect = 1,yreversed=true)
    ax_profile_x = GM.Axis(fig,title = "Ix")
    fig[2,:] = hgrid!(ax_Ex_xy,ax_Ey_xy,ax_Ez_xy,ax_Etot_xy,ax_profile_x)

    slider_ring1 = Slider(fig,color_active = :gray,halign = :right,range = 0:0.01:1, startvalue = 0,width = 500,linewidth = 20)
    slider_ring2 = Slider(fig,color_active = :gray,halign = :right,range = 0:0.01:1, startvalue = 0,width = 500,linewidth = 20)
    tb_ring1 = Textbox(fig,placeholder = "0.0", validator = Float64,height = 30,width = 40,fontsize = 12, tellwidth = false)
    tb_ring2 = Textbox(fig,placeholder = "0.0", validator = Float64,height = 30,width = 40,fontsize = 12, tellwidth = false)
    fig[3,:] = hgrid!(Label(fig,"ring1"),slider_ring1, tb_ring1)
    fig[4,:] = hgrid!(Label(fig,"ring2"),slider_ring2, tb_ring2)

    sz=40 
    cc = Int(sz/2)
    pos_emitter = (sz/2,sz/2,0.0)

    lift(slider_ring1.value,slider_ring2.value) do r1,r2
        tb_ring1.displayed_string = string(r1)
        tb_ring2.displayed_string = string(r2)
        fθ = (x) -> phase_theta(x, na, n[1], [r1,r2, 1.0])
        fϕ = (x) -> phase_phi(x)
        p=PSF.Focusfield(na,λ,n,pixelsize; ksize=128,excitationfield=[1.0,1*im],fθ=fθ, zstage=0.0) 

        roi=[(y,x,z) for y=1:sz, x=1:sz, z=0:0] 
        p.electricfield = 'x'
        imx_xy[] = PSF.pdf(p,roi,pos_emitter)[:,:,1]
        p.electricfield = 'y'
        imy_xy[] = PSF.pdf(p,roi,pos_emitter)[:,:,1]
        p.electricfield = 'z'
        imz_xy[] = PSF.pdf(p,roi,pos_emitter)[:,:,1]
        imall_xy[] = imx_xy[].+imy_xy[].+imz_xy[]
        I_x[] = imall_xy[][:,cc]
        

        roi=[(y,x,z) for y=sz/2:sz/2, x=1:sz, z=0:-0.04:-1.5] 

        p.electricfield = 'x'
        imx_xz[] = PSF.pdf(p,roi,pos_emitter)[1,:,:]
        p.electricfield = 'y'
        imy_xz[] = PSF.pdf(p,roi,pos_emitter)[1,:,:]
        p.electricfield = 'z'
        imz_xz[] = PSF.pdf(p,roi,pos_emitter)[1,:,:]
        imall_xz[] = imx_xz[].+imy_xz[].+imz_xz[]
        I_z[] = imall_xz[][cc,:]
    end

    GM.heatmap!(ax_Ex_xz,imx_xz,colormap = :inferno)
    GM.heatmap!(ax_Ey_xz,imy_xz,colormap = :inferno)
    GM.heatmap!(ax_Ez_xz,imz_xz,colormap = :inferno)
    GM.heatmap!(ax_Etot_xz,imall_xz,colormap = :inferno)
    GM.heatmap!(ax_Ex_xy,imx_xy,colormap = :inferno)
    GM.heatmap!(ax_Ey_xy,imy_xy,colormap = :inferno)
    GM.heatmap!(ax_Ez_xy,imz_xy,colormap = :inferno)
    GM.heatmap!(ax_Etot_xy,imall_xy,colormap = :inferno)
    GM.lines!(ax_profile_z,I_z)
    GM.lines!(ax_profile_x,I_x)

    GM.activate!(title="focus field")
    display(GM.Screen(), fig)

    return nothing
end

focus_gui(na,λ,n,pixelsize)


#look at pupil
hx = p.pupilfunctionx.pupil[:,:,1].*exp.(im*p.pupilfunctionx.pupil[:,:,2])
hy = p.pupilfunctiony.pupil[:,:,1].*exp.(im*p.pupilfunctiony.pupil[:,:,2])
hz = p.pupilfunctionz.pupil[:,:,1].*exp.(im*p.pupilfunctionz.pupil[:,:,2])
fig = Figure(resolution = (600, 400))
ax = CM.Axis(fig[1, 1],title="Ex",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, abs.(hx),colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 2],title="Ey",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, abs.(hy),colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 3],title="Ez",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, abs.(hz),colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 1],title="Ex",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, angle.(hx),colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 2],title="Ey",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, angle.(hy),colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 3],title="Ez",aspect=1,titlesize=20)
hm = CM.heatmap!(ax, angle.(hz),colormap=:inferno)
hidedecorations!(ax)
fig

#simulate PSF
ringsz = [0.0, 0.51, 1.0]
fθ = (x) -> phase_theta(x, na, n[1], ringsz) # 3D: [0.42, 0.61, 1.0]
fϕ = (x) -> phase_phi(x)
p=PSF.Focusfield(na,λ,n,pixelsize; ksize=128,excitationfield=[1.0,1*im],fθ=fθ,fϕ=fϕ,zstage=0.0) 
sz=50 
roi=[(y,x,z) for y=1:sz, x=1:sz, z=0:-0.04:-1.5] 

pos_emitter = (sz/2+1,sz/2+1,0.0)
p.electricfield = 'x'
imx=PSF.pdf(p,roi,pos_emitter)
p.electricfield = 'y'
imy=PSF.pdf(p,roi,pos_emitter)
p.electricfield = 'z'
imz=PSF.pdf(p,roi,pos_emitter)
imall = imx.+imy.+imz
cc = Int(pos_emitter[1])
#look at psf in x and y polarization
fig = Figure(resolution = (600, 400))
ax = CM.Axis(fig[1, 1],title="Ex psf",aspect=1)
hm = CM.heatmap!(ax, imx[cc,:,:],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 2],title="Ey psf",aspect=1)
hm = CM.heatmap!(ax, imy[cc,:,:],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 3],title="Ez psf",aspect=1)
hm = CM.heatmap!(ax, imz[cc,:,:],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 4],title="total psf",aspect=1)
hm = CM.heatmap!(ax, imall[cc,:,:],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 1],title="Ex psf",aspect=1)
hm = CM.heatmap!(ax, imx[:,:,1],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 2],title="Ey psf",aspect=1)
hm = CM.heatmap!(ax, imy[:,:,1],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 3],title="Ez psf",aspect=1)
hm = CM.heatmap!(ax, imz[:,:,1],colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[2, 4],title="total psf",aspect=1)
hm = CM.heatmap!(ax, imall[:,:,1],colormap=:inferno)
hidedecorations!(ax)
CM.activate!()
display(fig)

using HDF5

datadir = raw"C:\Users\Sheng\Documents\git\git_embl\MINFLUX-analysis\minflux demo\simdata/"
filename = "exfield_3d_642_2D_4.h5"
fh = h5open(datadir*filename, "w")
write(fh, "imx", imx)
write(fh, "imy", imy)
write(fh, "imz", imz)
write(fh, "imall", imall)
write(fh, "ringsz", ringsz)
write(fh, "phasemask", p.phasemask)
close(fh)


datadir = raw"C:\Users\Sheng\Documents\git\git_embl\MINFLUX-analysis\minflux demo\simdata/"
filename = "exfield_3d_642_2D_phasemask.h5"
fh = h5open(datadir*filename, "w")
write(fh, "phasemask", p.phasemask)
close(fh)


