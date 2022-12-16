# Create interpolated PSF from cubic spline

# Helper functions 


"""
    SplinePSF(psfstack) 

3D psf interpolated from a 3D PSF stack (spline interpolation)

#Fields
- `psfstack`             : a 3D PSF stack        
- `sitp`                 : spline interpolation object
- `pixelsize`            : Linear size of a back-projected pixel

"""
mutable struct SplinePSF{T<:AbstractFloat} <: PSF
    psfstack::Array
    sitp
    pixelsize::T

end

function SplinePSF(psfstack;pixelsize=0.1)

 
    itp = interpolate(psfstack, BSpline(Cubic(Line(OnGrid()))))
    psf_size = size(psfstack)
    x = -psf_size[1]/2+0.5:psf_size[1]/2-0.5
    y = -psf_size[2]/2+0.5:psf_size[2]/2-0.5
    z = -psf_size[3]/2+0.5:psf_size[3]/2-0.5
    sitp = scale(itp, x, y, z)


    return SplinePSF(psfstack,sitp,pixelsize)
end
# calculations
function pdf(p::SplinePSF, pixel::Tuple,x_emitter::Tuple)
    x=x_emitter.-pixel # convert to pixels
    return p.sitp[x...]
end    

   


# needed for broadcasting
Base.broadcastable(x::SplinePSF) = Ref(x)



