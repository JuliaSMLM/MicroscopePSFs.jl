
struct PupilFunction1d{T <: AbstractFloat} <: PSF 
    nₐ::T
    λ::T
    n::T
    pixelsize::T
    kpixelsize::T
    pupil::Array
    fintegeral::Function
end 

function pdfₐ(pupil::Array,kpixelsize,x,y,z,n,λ,f::Function)
        
    ksize = size(pupil, 2) # square 
    a_complex = ComplexF64(0.0)
    for ii in 1:ksize
        kr = kpixelsize * (ii - 0.5)
        kr2 = kr^2
        defocus = z * sqrt(ComplexF64((n / λ)^2 - kr2))
        a_complex += f(kr,x,y)*exp(2 * pi * im * defocus)*kr*kpixelsize #f(kr[ii],x,y)
    end
    #defocus = z .* sqrt.(complex((n / λ)^2 .- kr2))
    #a_complex = f.(kr,Ref(x),Ref(y)).*exp.(2 * pi * im .* defocus).*kr.*dkr
    #return a_complex*kpixelsize
    #return reduce(+,a_complex)
    return a_complex
end




function pdfₐ(p::PupilFunction1d, pixel::Tuple,x_emitter::Tuple)
    
    # calculate a defocus Phase
    #Δx=p.pixelsize.*(x_emitter[2].-pixel[2])
    #Δy=p.pixelsize.*(x_emitter[1].-pixel[1])
    Δx=p.pixelsize.*(pixel[2].-x_emitter[2])
    Δy=p.pixelsize.*(pixel[1].-x_emitter[1])

    Δz=x_emitter[3].-pixel[3]

    return pdfₐ(p.pupil,p.kpixelsize,Δx,Δy,Δz,p.n,p.λ,p.fintegeral)

end


function pdf(p::PupilFunction1d, pixel::Tuple,x_emitter::Tuple)
    out = pdfₐ(p::PupilFunction1d, pixel::Tuple,x_emitter::Tuple)
    return real(out*conj(out))
end
