# Create interpolated PSF from cubic spline

# Helper functions 


"""
    InterpolatedPSF(p,r;subsampling=4) 

3D psf interpolated from other PSF models 

#Fields
- psfstack             : a 3D psf stack        
- sitp                 : spline PSF object

"""
mutable struct SplinePSF <: PSF
    psfstack::Array
    sitp


end

function SplinePSF(psfstack)

    # setup interpolation 
    # we are working in a 'pixel' 
    # unit basis 



    #integration over finite pixels
    #TODO!


 
    itp = interpolate(psfstack, BSpline(Cubic(Line(OnGrid()))))
    psf_size = size(psfstack)
    x = -psf_size[1]/2+0.5:psf_size[1]/2-0.5
    y = -psf_size[2]/2+0.5:psf_size[2]/2-0.5
    z = -psf_size[3]/2+0.5:psf_size[3]/2-0.5
    sitp = scale(itp, x, y, z)


    return SplinePSF(psfstack,sitp)
end
# calculations
function pdf(p::SplinePSF, pixel::Tuple,x_emitter::Tuple)
    x=x_emitter.-pixel # convert to pixels
    return p.sitp[x...]
end    

   


# needed for broadcasting
Base.broadcastable(x::SplinePSF) = Ref(x)



