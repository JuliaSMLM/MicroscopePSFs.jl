## Scalar approximation of high NA objective PSFs 

# Helper functions 


"""
    ImmPSF 

3D psf with index mismatch aberration using scalar or vector model and OTF rescaling 

# Fields
- `pupilfunction`   : Pupil Function structure        
- `pixelsize`       : Linear size of a back-projected pixel, unit: micron  
- 'Σ'               : OTF rescaling via image space 2D Gaussian Covariance matrix 
- 'ksize'           : number of pixels in pupil 

# Constructor
    ImmPSF(nₐ, λ, n, pixelsize; inputpupil=nothing, Σ=0, ksize=256, z::ZernikeCoefficients=ZernikeCoefficients(1))
- 'nₐ'          : numerical aperture
- 'λ'           : emission wavelength, unit is micrion
- 'n'           : refractive indices of (sample medium, coverglass, immersion medium)
- 'zstage'      : position of the sample stage, unit: micron
- 'inputpupil'  : input pupil function, default is nothing
- 'mvtype'      : type of motion, the z position is changed by moving 'bead' or 'stage'
"""
mutable struct ImmPSF{PF<:PupilFunction,T<:AbstractFloat,I<:Int} <: PSF
    pupilfunction::Vector{PF}
    pixelsize::T
    Σ::T
    ksize::Int
end

function ImmPSF(nₐ, λ, n::Vector, pixelsize; zstage = 0.0,inputpupil=nothing, Σ=0, ksize=256, z::ZernikeCoefficients=ZernikeCoefficients(1),mvtype="bead")

    pupil = [zeros(ksize, ksize, 2) for x=1:6]
    kpixelsize = 2 * nₐ / λ / ksize
    
    k0 = (ksize + 1) / 2

    for ii in 1:ksize, jj in 1:ksize # jj is index along y  
        kx = kpixelsize * (ii - k0)
        ky = kpixelsize * (jj - k0)
        
        kr2 = kx^2 + ky^2
        Tp, Ts, sinθ₁, cosθ₁, _, cosθ₃ = calFresnel(kr2,λ,n)
        #immphase = exp(2*pi*(n[1]/λ*cosθ₁*n[1]/n[3]*zstage-n[3]/λ*cosθ₃*zstage)*im)
        immphase = exp(-2*pi*(n[3]/λ*cosθ₃*zstage)*im)

        if kr2 < (nₐ / λ)^2 
            ρ=sqrt(kr2)/(nₐ / λ)
            ϕ=atan(ky,kx)
            _, _, h = calEfield(ϕ, Tp, Ts, sinθ₁, cosθ₁)
            pupil_mag = 0.0
            pupil_phase = 0.0
            for nn=1:length(z.mag)
                if abs(z.mag[nn])>0.0
                    pupil_mag += z.mag[nn]*zernikepolynomial(nn-1,ρ,ϕ)
                end
            end
            pupil_mag *= abs(immphase)

            for nn=1:length(z.phase)
                if abs(z.phase[nn])>0.0
                    pupil_phase += z.phase[nn]*zernikepolynomial(nn-1,ρ,ϕ)
                end
            end
            pupil_phase += angle(immphase)
            for l in eachindex(h)
                pupil[l][jj,ii,1] = pupil_mag*abs(h[l])
                pupil[l][jj,ii,2] = pupil_phase+angle(h[l])
            end

        end
    end


    normf = 0.0
    for j=eachindex(pupil)
        normf += sum(pupil[j][:,:,1].^2)
        if inputpupil !== nothing
            pupil[j][:,:,1] .*= inputpupil[:,:,1]
            pupil[j][:,:,2] .+= inputpupil[:,:,2]
        end
        #normf += sum(pupil[j][:,:,1].^2)
    end

    normf = sqrt(normf)/kpixelsize
    for j=eachindex(pupil)
        pupil[j][:,:,1] = pupil[j][:,:,1]./normf
    end
    

    p=Vector{PupilFunction}()
    for j=eachindex(pupil)
        if mvtype == "bead"
            append!(p,[PupilFunction(nₐ, λ, n[1], pixelsize, kpixelsize, pupil[j])])
        end
        if mvtype == "stage"
            append!(p,[PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupil[j])])
        end
    end
   
    return ImmPSF{PupilFunction,typeof(nₐ),Int}(p, pixelsize, Σ, ksize)
end

function pdfₐ(p::ImmPSF, pixel::Tuple,x_emitter::Tuple)

        out=pdfₐ.(p.pupilfunction, Ref(pixel), Ref(x_emitter)) * p.pixelsize
    
    return out
end    

function pdfₐ(p::ImmPSF, roi::Array,x_emitter::Tuple)
    out=[zeros(Complex,size(roi)) for x=1:length(p.pupilfunction)]
    Threads.@threads for ii=1:length(roi)
        outpx=pdfₐ(p, roi[ii], x_emitter)
        for j=eachindex(outpx)
            out[j][ii] = outpx[j]
        end
    end
    return out
end    

function pdf(p::ImmPSF, pixel::Tuple,x_emitter::Tuple)
    return sum(abs2.(pdfₐ(p, pixel, x_emitter))) 
end    

function pdf(p::ImmPSF,roi::Array,x_emitter::Tuple)
    out=zeros(size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdf(p,roi[ii],x_emitter)
    end
    return out
end


