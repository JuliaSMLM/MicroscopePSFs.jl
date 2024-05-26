# Create interpolated PSF from cubic spline

# Helper functions 


"""
    SplinePSF

3D psf interpolated from a 3D PSF stack (spline interpolation)

# Fields
- `psfstack`             : a 3D PSF stack        
- `sitp`                 : spline interpolation object
- `pixelsize`            : Linear size of a back-projected pixel
- `pixelsize_z`          : step size in z, unit: micron

# Constructor
    SplinePSF(psfstack;pixelsize_z=0.05,pixelsize=0.1)
"""
mutable struct SplinePSF{T<:AbstractFloat} <: PSF
    psfstack::Array
    sitp
    pixelsize::T
    pixelsize_z::T
end

function SplinePSF(psfstack::Array{Float32,3};pixelsize_z=0.05,pixelsize=0.1)

 
    itp = interpolate(psfstack, BSpline(Cubic(Line(OnGrid()))))
    psf_size = size(psfstack)
    x = -psf_size[1]/2+0.5:psf_size[1]/2-0.5
    y = -psf_size[2]/2+0.5:psf_size[2]/2-0.5
    z = (-psf_size[3]/2+0.5:psf_size[3]/2-0.5)*pixelsize_z
    sitp = scale(itp, x, y, z)

    return SplinePSF(psfstack,sitp,pixelsize,pixelsize_z)
end

function SplinePSF(psfstack::Array{Float32,4};pixelsize_z=0.05,pixelsize=0.1,pixelsize_st=0.01)

 
    itp = interpolate(psfstack, BSpline(Cubic(Line(OnGrid()))))
    psf_size = size(psfstack)
    y = -psf_size[1]/2+0.5:psf_size[1]/2-0.5
    x = -psf_size[2]/2+0.5:psf_size[2]/2-0.5
    z = (0:psf_size[3]-1)*pixelsize_z
    zs = (0:psf_size[4]-1)*pixelsize_st
    sitp = scale(itp, y, x, z, zs)

    return SplinePSF(psfstack,sitp,pixelsize,pixelsize_z)
end

function SplinePSF(psfstack::Array{Float32,5};pixelsize_z=0.05,pixelsize=0.1,image_size=[256,256])

 
    itp = interpolate(psfstack, BSpline(Cubic(Line(OnGrid()))))
    psf_size = size(psfstack)
    x = -psf_size[1]/2+0.5:psf_size[1]/2-0.5
    y = -psf_size[2]/2+0.5:psf_size[2]/2-0.5
    z = (0:psf_size[3]-1)*pixelsize_z
    corx = range(0,image_size[1],psf_size[4])
    cory = range(0,image_size[2],psf_size[5])
    sitp = scale(itp, x, y, z, corx, cory)
   
    return SplinePSF(psfstack,sitp,pixelsize,pixelsize_z)
end
# calculations
function pdf(p::SplinePSF, pixel::Tuple,x_emitter::Tuple)
    x=x_emitter.-pixel # convert to pixels
    return p.sitp[x...]
end    

function pdf(p::SplinePSF, pixel::Tuple,x_emitter::Tuple, cor::Tuple)  
    cor_emitter = cor.+x_emitter[1:2]
    x=x_emitter.-pixel # convert to pixels
    x = (x...,cor_emitter...)
    return p.sitp[x...]
end   

function pdf(p::SplinePSF, pixel::Tuple,x_emitter::Tuple,z_stage::Float64)
    x=pixel.-x_emitter # convert to pixels
    x = (x[1],x[2],-x[3],z_stage)
    return p.sitp[x...]
end 
# needed for broadcasting
Base.broadcastable(x::SplinePSF) = Ref(x)



