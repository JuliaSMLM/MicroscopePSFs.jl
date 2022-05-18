## Abstract type and methods assuming a known pupil function 

"""
    PupilFunction

A structure used by 3D PSF methods.

#Fields
- `nₐ`              : Numerical Aperture           
- `λ`               : Emission Wavelength           
- `n`               : Immersion Index           
- `pixelsize`       : Linear size of a back-projected pixel
- 'kpixelsize'      : Number of pixels used in pupil function image
- 'pupil'           : Pupil function [ksize*ksize*2] (Amplitude, Phase)  

Calculations of various 3D PSF types will involve calculations based around 
Pupil Funciton.  They differ in the complexity and approach used to calculate 
Pupil Function. 
"""
mutable struct PupilFunction{T <: AbstractFloat} <: PSF 
    nₐ::T
    λ::T
    n::T
    pixelsize::T
    kpixelsize::T
    pupil::Array{T}
end 

function pdfₐ(pupil::Array{<:Real,3},kpixelsize,x,y,z,n,λ,nₐ)
    
    a_real = 0
    a_im = 0
    
    sz = size(pupil, 1) # square 

    # calculate a defocus Phase
    # A(x,y,z)=∱∱ pupil*exp(2πi(kx*x+ky*y))*exp( 2πi*z* sqrt( (n/λ)²-(kx²+ky²) ) )dkx dky
    
    kmax2 = (nₐ / λ)^2
    k0 = (sz + 1) / 2

    for ii in 1:sz, jj in 1:sz # jj is index along y  
        kx = kpixelsize * (ii - k0)
        ky = kpixelsize * (jj - k0)
        kr2 = kx^2 + ky^2
        
        if kr2 < kmax2
            
            idx1 = (ii - 1) * sz + jj
            idx2 = sz * sz + (ii - 1) * sz + jj

            defocus = z * sqrt((n / λ)^2 - kr2)
            θ = pupil[idx2] + 2 * pi * (defocus + x * kx + y * ky)
            a_real += pupil[idx1] * cos(θ)
            a_im += pupil[idx1] * sin(θ)
        end

    end

    return complex(a_real, a_im)

end



function pdfₐ(p::PupilFunction, pixel::Tuple,x_emitter::Tuple)
    
    # calculate a defocus Phase
    ΔX=p.pixelsize.*(x_emitter.-pixel)
    
    return pdfₐ(p.pupil,p.kpixelsize,ΔX[1],ΔX[2],ΔX[3],p.n,p.λ,p.nₐ)

end


function pdf(p::PupilFunction, pixel::Tuple,x_emitter::Tuple)
    return abs2(pdfₐ(p::PupilFunction, pixel::Tuple,x_emitter::Tuple))
end



"""
    normalize!(p::PupilFunction)

normalize pupil using Parseval's theorem
"""
function normalize!(p::PupilFunction)

    sz = size(p.pupil, 1) # square 
    α = 0

    # find total energy
    for ii in 1:sz, jj in 1:sz # jj is index along y  
        idx1 = (ii - 1) * sz + jj
        α += p.pupil[idx1]^2
    end

    α = sqrt(α)/p.kpixelsize

     # normalize
    for ii in 1:sz, jj in 1:sz # jj is index along y  
        idx1 = (ii - 1) * sz + jj
        p.pupil[idx1] = p.pupil[idx1] / α
    end

end


