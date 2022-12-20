## Scalar approximation of high NA objective PSFs 

# Helper functions 


"""
    Scalar3D 

3D psf using scalar model and OTF rescaling 

# Fields
- `pupilfunction`   : Pupil Function structure           
- `pixelsize`            : Linear size of a back-projected pixel
- 'Σ'               : OTF rescaling via image space 2D Gaussian Covariance matrix 
- 'ksize'           : number of pixels in pupil 


"""
mutable struct Scalar3D{PF<:PupilFunction,T<:AbstractFloat,I<:Int} <: PSF
    pupilfunction::PF
    pixelsize::T
    Σ::T
    ksize::Int
end

function Scalar3D(nₐ, λ, n, pixelsize;inputpupil=nothing, Σ=0, ksize=256, z::ZernikeCoefficients=ZernikeCoefficients(1))

    pupil = zeros(ksize, ksize, 2)
    kpixelsize = 2 * nₐ / λ / ksize
    
    k0 = (ksize + 1) / 2

    for ii in 1:ksize, jj in 1:ksize # jj is index along y  
        kx = kpixelsize * (ii - k0)
        ky = kpixelsize * (jj - k0)
        
        kr2 = kx^2 + ky^2
        if kr2 < (nₐ / λ)^2 
            ρ=sqrt(kr2)/(nₐ / λ)
            ϕ=atan(ky,kx)
            for nn=1:length(z.mag)
                if abs(z.mag[nn])>0.0
                    pupil[jj,ii,1] += z.mag[nn]*zernikepolynomial(nn-1,ρ,ϕ)
                end
            end

            for nn=1:length(z.phase)
                if abs(z.phase[nn])>0.0
                    pupil[jj,ii,2] += z.phase[nn]*zernikepolynomial(nn-1,ρ,ϕ)
                end
            end

        end
    end
    if inputpupil !== nothing
        pupil[:,:,1] .*= inputpupil[:,:,1]
        pupil[:,:,2] .+= inputpupil[:,:,2]
    end
    p = PupilFunction(nₐ, λ, n, pixelsize, kpixelsize, pupil)
    normalize!(p)
    return Scalar3D{PupilFunction,typeof(nₐ),Int}(p, pixelsize, Σ, ksize)
end

function pdfₐ(p::Scalar3D, pixel::Tuple,x_emitter::Tuple)
    return pdfₐ(p.pupilfunction, pixel, x_emitter) * p.pixelsize 
end    

function pdfₐ(p::Scalar3D, roi::Array,x_emitter::Tuple)
    out=zeros(Complex,size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdfₐ(p, roi[ii], x_emitter)
    end
    return out
end    

function pdf(p::Scalar3D, pixel::Tuple,x_emitter::Tuple)
    return abs2(pdfₐ(p, pixel, x_emitter)) 
end    

function pdf(p::Scalar3D,roi::Array,x_emitter::Tuple)
    out=zeros(size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdf(p,roi[ii],x_emitter)
    end
    return out
end


