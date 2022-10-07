## Scalar approximation of high NA objective PSFs 

# Helper functions 


"""
    Dipole3D 

3D psf using vector model and OTF rescaling 

#Fields
- `pupilfunction`   : Pupil Function structure           
- 'Σ'               : OTF rescaling via image space 2D Gaussian Covariance matrix 
- 'ksize'           : number of pixels in pupil 
- 'dipole_ang'         : orientation of dipole moment (polar, azimuthal) 

"""
mutable struct Dipole3D{PF<:PupilFunction,T<:AbstractFloat,I<:Int} <: PSF
    pupilfunctionx::PF
    pupilfunctiony::PF
    pixelsize::T
    Σ::T
    dipole_ang::Vector
    ksize::Int
    electricfield::Char
end

"""
    Dipole3D()
    - 'n'       : refractive indices (sample medium, cover glass, immersion medium)
"""
function Dipole3D(nₐ, λ, n::Vector, pixelsize, dipole_ang::Vector; electricfield='x', Σ=0.0, ksize=256, z::ZernikeCoefficients=ZernikeCoefficients(1))

    pupilx = zeros(ksize, ksize, 2)
    pupily = zeros(ksize, ksize, 2)
    kpixelsize = 2 * nₐ / λ / ksize
    
    k0 = (ksize + 1) / 2

    α = dipole_ang[1]
    β = dipole_ang[2]
    dvec = [sin(α)*cos(β),sin(α)*sin(β),cos(α)]
    for ii in 1:ksize, jj in 1:ksize # jj is index along y  
        kx = kpixelsize * (ii - k0)
        ky = kpixelsize * (jj - k0)
        
        kr2 = kx^2 + ky^2
    

        Tp, Ts, sinθ₁, cosθ₁, _, _=calFresnel(kr2,λ,n)
        
        if kr2 < (nₐ / λ)^2 
            ρ=sqrt(kr2)/(nₐ / λ)
            ϕ=atan(ky,kx)

            hx, hy, _ = calEfield(ϕ, Tp, Ts, sinθ₁, cosθ₁; dvec=dvec)
            for nn=1:length(z.mag)
                if abs(z.mag[nn])>0.0
                    pupilx[jj,ii,1] += z.mag[nn]*zernikepolynomial(nn-1,ρ,ϕ)
                    
                end
            end
            pupily[jj,ii,1] = pupilx[jj,ii,1]
            pupilx[jj,ii,1] *= abs(hx)
            pupily[jj,ii,1] *= abs(hy)
            for nn=1:length(z.phase)
                if abs(z.phase[nn])>0.0
                    pupilx[jj,ii,2] += z.phase[nn]*zernikepolynomial(nn-1,ρ,ϕ)
                end
            end
            pupily[jj,ii,2] = pupilx[jj,ii,2]
            pupilx[jj,ii,2] += angle(hx)
            pupily[jj,ii,2] += angle(hy)
        end
    end

    #pupilx = cat(abs.(Ex),angle.(Ex).+pupilphase;dims=3)
    #pupily = cat(abs.(Ey),angle.(Ey).+pupilphase;dims=3)

    normf = sqrt(sum(pupilx[:,:,1].^2) + sum(pupily[:,:,1].^2))/kpixelsize


    pupilx[:,:,1] = pupilx[:,:,1]./normf
    pupily[:,:,1] = pupily[:,:,1]./normf

    px = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupilx)
    py = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupily)


    #normalize!(px)
    #normalize!(py)
    return Dipole3D{PupilFunction,typeof(nₐ),Int}(px,py, pixelsize, Σ,dipole_ang,ksize,electricfield)
end

function pdfₐ(p::Dipole3D, pixel::Tuple,x_emitter::Tuple)
    if p.electricfield == 'x'
        out = pdfₐ(p.pupilfunctionx, pixel, x_emitter) * p.pixelsize 
    end
    if p.electricfield == 'y'
        out = pdfₐ(p.pupilfunctiony, pixel, x_emitter) * p.pixelsize 
    end
    return out
end    

function pdfₐ(p::Dipole3D, roi::Array,x_emitter::Tuple)
    out=zeros(Complex,size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdfₐ(p, roi[ii], x_emitter)
    end
    return out
end    

function pdf(p::Dipole3D, pixel::Tuple,x_emitter::Tuple)
    return abs2(pdfₐ(p, pixel, x_emitter)) 
end    

function pdf(p::Dipole3D,roi::Array,x_emitter::Tuple)
    out=zeros(size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdf(p,roi[ii],x_emitter)
    end
    return out
end


function calFresnel(kr2,λ,n::Vector)
    sinθ₁ = sqrt(kr2)*λ/n[1]
    cosθ₁ = sqrt(complex(1-kr2*λ*λ/n[1]/n[1])) 
    cosθ₂ = sqrt(complex(1-kr2*λ*λ/n[2]/n[2]))
    cosθ₃ = sqrt(complex(1-kr2*λ*λ/n[3]/n[3]))

    FresnelP₁₂ = 2*n[1]*cosθ₁/(n[1]*cosθ₂+n[2]*cosθ₁)
    FresnelS₁₂ = 2*n[1]*cosθ₁/(n[1]*cosθ₁+n[2]*cosθ₂)  
    FresnelP₂₃ = 2*n[2]*cosθ₂/(n[2]*cosθ₃+n[3]*cosθ₂) 
    FresnelS₂₃ = 2*n[2]*cosθ₂/(n[2]*cosθ₂+n[3]*cosθ₃) 
    Tp = FresnelP₁₂*FresnelP₂₃
    Ts = FresnelS₁₂*FresnelS₂₃

    return Tp, Ts, sinθ₁, cosθ₁, cosθ₂, cosθ₃
end


function calEfield(ϕ, Tp, Ts, sinθ₁, cosθ₁; dvec = [1,1,1])

    pvec = [cosθ₁*cos(ϕ),cosθ₁*sin(ϕ),-sinθ₁].*Tp.*dvec
    svec = [-sin(ϕ),cos(ϕ),0.0].*Ts.*dvec
    hx = pvec.*cos(ϕ)-svec.*sin(ϕ)
    hy = pvec.*sin(ϕ)+svec.*cos(ϕ)
    return sum(hx), sum(hy), hcat(hx,hy)
end
