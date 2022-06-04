## Types and functions for 2D Gaussian psf


"""
    Gauss2D

Isotropic 2D Gaussian psf

#Fields
- `σ`            : Gaussian σ           
- `pixelsize`    : Linear size of a pixel

"""
mutable struct Gauss2D{T<:AbstractFloat} <: PSF
    σ::T 
    pixelsize::T
end

function pdfₐ(p::Gauss2D,pixel::Tuple,x_emitter::Tuple)
    return sqrt(pdf(p,pixel,x_emitter))  
end    

function pdf(p::Gauss2D,pixel::Tuple,x_emitter::Tuple)
    r²=((x_emitter[1]-pixel[1]).^2+
    (x_emitter[2]-pixel[2]).^2)*p.pixelsize^2
   
    return p.pixelsize^2/(2*pi*p.σ^2)*
        exp(-r²/(2*p.σ^2))  
end    

#conversion functions
"""
    Gauss2D(p::Airy2D)

convert an Airy2D to Gaussian PSF    
"""
function Gauss2D(p::Airy2D)
    ν=2*pi*p.nₐ/p.λ
    σ = 0.42*π/ν
    return Gauss2D(σ,p.pixelsize)
end





