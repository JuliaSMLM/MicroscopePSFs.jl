# Create interpolated PSF

# Helper functions 


"""
    InterpolatedPSF(p,r;subsampling=4) 

3D psf interpolated from other PSF models 

#Fields
- p             : Generator PSF           
- r             : range (x,y,z) where x,y,z are maximum supported range in each dimension
- subsampling   : subsampling (default = 4)

"""
mutable struct InterpolatedPSF <: PSF
    p::PSF 
    r
    subsampling   
    itp
    itp_real
    itp_imag
    
function InterpolatedPSF(p,r; subsampling=4)

    # setup interpolation 
    # we are working in a 'pixel' 
    # unit basis 

    n=length(r)
    #sampling the PSF
    dx=1/subsampling
    dz=p.pixelsize/subsampling #could make better 

    x=Array(-r[1]:dx:r[1])
    y=Array(-r[2]:dx:r[2])
    x_emitter=(0,0)
    if n==2
        ndrange=(x,y)
        roi=[(i,j) for i=x, j=y]
    end

    if n==3
        z=Array(-r[3]:dz:r[3])
        ndrange=(x,y,z)
        roi=[(i,j,k) for j=y, i=x, k=z]
        x_emitter=(0,0,0)
    end

    #integration over finite pixels
    #TODO!

    im=pdf(p,roi,x_emitter)
    itp=interpolate(ndrange,im,Gridded(Linear()))

    im=pdfₐ(p,roi,x_emitter)
    itp_real=interpolate(ndrange,real.(im),Gridded(Linear()))
    itp_imag=interpolate(ndrange,imag.(im),Gridded(Linear()))

    return new(p, r, subsampling, itp,itp_real,itp_imag)
end

end

# calculations
function pdf(p::InterpolatedPSF, pixel::Tuple,x_emitter::Tuple)
    x=x_emitter.-pixel # convert to pixels
    return p.itp[x...]
end    

function pdfₐ(p::InterpolatedPSF, pixel::Tuple,x_emitter::Tuple)
    x=x_emitter.-pixel # convert to pixels
    
    return Complex(p.itp_real[x...],p.itp_imag[x...])
end    


# needed for broadcasting
Base.broadcastable(x::InterpolatedPSF) = Ref(x)



