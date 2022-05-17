## Types and functions for 2D Gaussian psf



"""
    Airy PSF 

2D Airy Patter psf using paraxial, scalar model

#Fields
- `na`              : Numerical Aperture           
- `λ`               : Emission Wavelength           
- `pixelsize`       : Linear size of a back-projected pixel

The Airy PSF is  
I(r)=ν²/(4π)(2*J₁(ν*r)/(ν*r))²    
where     
ν=πD/(λf)=2*π*nₐ/λ  


"""
mutable struct Airy2D{T<:AbstractFloat} <: PSF
    nₐ::T
    λ::T 
    pixelsize::T
    ν::T
    function Airy2D(nₐ,λ,pixelsize)
        ν=2*pi*nₐ/λ
        return new{typeof(nₐ)}(nₐ,λ,pixelsize,ν)    
    end
end


function pdfₐ(p::Airy2D,x_pixel::Tuple,x_emitter::Tuple)
    r=sqrt((x_emitter[1]-x_pixel[2]).^2+
    (x_emitter[2]-x_pixel[1]).^2)*p.pixelsize
    w = r * p.ν
    w = max(w, 1f-5) #handle zero case   
    return p.pixelsize*p.ν/(sqrt(4π))*(2*besselj1(w)/(w))
end    

function pdf(p::Airy2D,x_pixel::Tuple,x_emitter::Tuple)
    return abs2(pdfₐ(p,x_pixel,x_emitter))
end    









