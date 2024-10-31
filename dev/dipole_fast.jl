# 1D integral

using Revise
using MicroscopePSFs
PSF=MicroscopePSFs
using CairoMakie
CM = CairoMakie
using BenchmarkTools
using SpecialFunctions
# Create a scalar PSF
na=1.5
n=[1.33,1.516,1.516] # refractive indices (sample medium, cover glass, immersion)
λ=.69
pixelsize=λ/na/16
dipole_ang = [90,0].*pi./180
p=PSF.Dipole3D_fast(na,λ,n,pixelsize,dipole_ang; ksize=256,excitationfield=[1.0,0.0,0.0],mvtype="stage",δ=nothing); # scalar excitation, for fast rotating dipole

#simulate PSF
sz=200
roi=[(y,x,0) for y=1:sz, x=1:sz] 

pos_emitter = (sz/2,sz/2,0.0)
p.electricfield = 'x'
@time imx=PSF.pdf(p,roi,pos_emitter);
@profview imx=PSF.pdf(p,roi,pos_emitter);
p.electricfield = 'y'
@time imy=PSF.pdf(p,roi,pos_emitter);
@profview imy=PSF.pdf(p,roi,pos_emitter);


#look at psf in x and y polarization
fig = Figure(size = (800, 400))
ax = CM.Axis(fig[1, 1],title="Ex psf",aspect=1)
hm = CM.heatmap!(ax, imx,colormap=:inferno)
hidedecorations!(ax)
ax = CM.Axis(fig[1, 2],title="Ey psf",aspect=1)
hm = CM.heatmap!(ax, imy,colormap=:inferno)
hidedecorations!(ax)
fig

# kmax = na / λ
# kr = (range(0,1,length=p.ksize+1)).*kmax
# dkr = diff(kr)
# kr = kr[1:end-1].+dkr./2
# a_complex = zeros(ComplexF64, p.ksize)
# for ii in eachindex(kr)
#     a_complex[ii] = p.pupilfunctionx.fintegeral(kr[ii],0.0,0.0) *kr[ii]*dkr[ii]
# end

# sum(a_complex)


f = (x) -> PSF.pdf(p,roi[1],(x[1],x[2],x[3]))*x[4]+x[5] 

using Zygote
#grad_s = Zygote.gradient((x) -> f(x), 10.5) 

@time grad_s = Zygote.forward_jacobian((x) -> f(x), [10.5,10.5,0.0,1.0,0.0])
# using CUDA
#     # Convert inputs to CuArrays
#     sz = 200*16
# out_d = CUDA.zeros(sz*sz)




# function fintegral(kr,x,y)
#     #λ=.69
#     #n = [1.33,1.52,1.52]
#     #kr2=kr^2
#     #dvec = [1.0,0.0,0.0]
#     #Tp, Ts, sinθ₁, cosθ₁,  _, cosθ₃ = PSF.calFresnel(kr2,λ,n)
#     #cos2ϕᵣ,sin2ϕᵣ,cosϕᵣ, sinϕᵣ, r = PSF.car2pol(x,y)
#     r = sqrt(x^2 + y^2)
#     sinϕᵣ = y/r
#     cosϕᵣ = x/r
#     cos2ϕᵣ = 2*cosϕᵣ^2-1
#     sin2ϕᵣ = 2*sinϕᵣ*cosϕᵣ

#     g = kr*r
#     J0 = SpecialFunctions.besselj(0,g)
#     J1 = SpecialFunctions.besselj(1,g)
#     J2 = SpecialFunctions.besselj(2,g)

#     Tp = 1.0
#     Ts = 1.0
#     cosθ₁ = 1.0-kr^2
#     sinθ₁ = kr
#     hx = (cosθ₁*pi*(J0-J2*cos2ϕᵣ)-cosθ₁*pi*J2*sin2ϕᵣ-sinθ₁*2*pi*J1*cosϕᵣ)*Tp + (pi*(J0+J2*cos2ϕᵣ)+pi*J2*sin2ϕᵣ)*Ts

#     #hx = pi*(J0-J2*cos2ϕᵣ)-pi*J2*sin2ϕᵣ-2*pi*J1*cosϕᵣ+pi*(J0+J2*cos2ϕᵣ)+pi*J2*sin2ϕᵣ
    
#     return hx
# end

# # Define the kernel function
# function kernel_pdfₐ(out)
#     na = 1.5
#     λ = .69
#     y = threadIdx().x
#     x = blockIdx().x
#     ii = y + (x - 1) * blockDim().x
#     kr = na/λ
#     #out[ii] = SpecialFunctions.besselj(0,r)
#     out[ii]= fintegral(kr,x,y)
#     return
# end

# # Launch the kernel
# threads_per_block = 1024
# blocks_per_grid = sz*sz ÷ threads_per_block
# @time @cuda threads=threads_per_block blocks=blocks_per_grid kernel_pdfₐ(out_d)

# # Copy the result back to the host
# out = Array(out_d)
# out = reshape(out, threads_per_block, blocks_per_grid)

# sz = 200*16
# out = zeros(sz,sz)
# roi=[(y,x) for y=1:sz, x=1:sz] 
# @time Threads.@threads for j in eachindex(roi)
#     x = roi[j][1]
#     y = roi[j][2]
#     r = (roi[j][1]^2 + roi[j][2]^2)/1000.0
#     out[j]= real(p.pupilfunctionx.fintegeral(r,x,y))
#     #out[j] = SpecialFunctions.besselj(0,r)
# end

xs = collect(1:sz/2)
fs = zeros(length(xs))
for j in eachindex(xs)

    x = xs[j]*pixelsize
    y = x
    r = sqrt(x^2 + y^2)
    cosϕᵣ = x/r
    a_complex = ComplexF64(0.0)
    ksize = size(p.pupilfunctionx.pupil, 2) 
    for ii in 1:ksize
        kr = p.pupilfunctionx.kpixelsize * (ii - 0.5)
        sinθ₁ = kr*λ/n[1]
        J1 = PSF.bJ(1,kr*r)
        a_complex += sinθ₁*2*pi*im*J1*cosϕᵣ*p.pupilfunctionx.kpixelsize
    end
    fs[j] = a_complex*conj(a_complex)
end