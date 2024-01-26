## Scalar approximation of high NA objective PSFs 

# Helper functions 


"""
    Dipole3D 

3D psf using vector model and OTF rescaling 

# Fields
- `pupilfunctionx`   : Pupil Function defines the x component of the electric field    
- `pupilfunctiony`   : Pupil Function defines the y component of the electric field
- `pixelsize`        : Linear size of a back-projected pixel, unit: micron       
- 'Σ'                : OTF rescaling via image space 2D Gaussian Covariance matrix 
- 'dipole_ang'       : orientation of dipole moment (polar, azimuthal) 
- 'ksize'            : number of pixels in pupil 
- 'electricfield'    : electric field component, 'x' or 'y'
- 'normf'            : normalization factor

# Constructor
    Dipole3D(nₐ, λ, n, pixelsize, dipole_ang; Σ=0, ksize=256, z::ZernikeCoefficients=ZernikeCoefficients(1))
- 'nₐ'          : numerical aperture
- 'λ'           : emission wavelength, unit is micrion
- 'n'           : refractive indices (sample medium, cover glass, immersion medium)
"""
mutable struct Dipole3D{PF<:PupilFunction,T<:AbstractFloat,I<:Int} <: PSF
    pupilfunctionx::PF
    pupilfunctiony::PF
    apodization::AbstractArray
    pixelsize::T
    Σ::T
    dipole_ang::Vector
    ksize::Int
    electricfield::Char
    normf::T
end


function Dipole3D(nₐ, λ, n::Vector, pixelsize, dipole_ang::Vector; electricfield='x', Σ=0.0, ksize=256, z::ZernikeCoefficients=ZernikeCoefficients(1))

    pupilx = zeros(ksize, ksize, 2)
    pupily = zeros(ksize, ksize, 2)
    apod = zeros(ksize,ksize,2)
    kpixelsize = 2 * nₐ / λ / ksize
    
    k0 = (ksize + 1) / 2

    α = dipole_ang[1]
    β = dipole_ang[2]
    dvec = [sin(α)*cos(β),sin(α)*sin(β),cos(α)]
    for ii in 1:ksize, jj in 1:ksize # jj is index along y  
        kx = kpixelsize * (ii - k0)
        ky = kpixelsize * (jj - k0)
        
        kr2 = kx^2 + ky^2
    

        Tp, Ts, sinθ₁, cosθ₁,  cosθ₂, cosθ₃ =calFresnel(kr2,λ,n)

        if kr2 < (nₐ / λ)^2 
            apod[jj,ii,1] = abs(sqrt(cosθ₃)/cosθ₁)
            apod[jj,ii,2] = angle(sqrt(cosθ₃)/cosθ₁)
   
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

    #normf = sqrt(sum(pupilx[:,:,1].^2) + sum(pupily[:,:,1].^2))/kpixelsize
    normf = pi/4/kpixelsize*ksize

    pupilx[:,:,1] = pupilx[:,:,1]./normf
    pupily[:,:,1] = pupily[:,:,1]./normf

    px = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupilx)
    py = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupily)

    normf1 = normf*kpixelsize/ksize
    #normalize!(px)
    #normalize!(py)
    return Dipole3D{PupilFunction,typeof(nₐ),Int}(px,py,apod, pixelsize, Σ,dipole_ang,ksize,electricfield,normf1)
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

"""
    calFresnel(kr2,λ,n::Vector)

Calculate Fresnel coefficients

# Arguments
- `kr2`     : magnitude square of the radial component of the k vector
- `λ`       : emission wavelength, unit is micrion
- `n`       : refractive indices (sample medium, cover glass, immersion medium)

# Returns
- `Tp`      : Fresnel coefficient for p-polarization
- `Ts`      : Fresnel coefficient for s-polarization
- `sinθ₁`   : sine of the incident angle in the sample medium
- `cosθ₁`   : cosine of the incident angle in the sample medium
- `cosθ₂`   : cosine of the incident angle in the cover glass
- `cosθ₃`   : cosine of the incident angle in the immersion medium
"""
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

"""
    calEfield(ϕ, Tp, Ts, sinθ₁, cosθ₁; dvec = [1,1,1])

Calculate electric field of a given dipole orientation

# Arguments
- `ϕ`       : azimuthal coordinate of the electric field
- `Tp`      : Fresnel coefficient for p-polarization
- `Ts`      : Fresnel coefficient for s-polarization
- `sinθ₁`   : sine of the incident angle in the sample medium
- `cosθ₁`   : cosine of the incident angle in the sample medium
- `dvec`    : dipole moment vector

# Returns
- x component of the electric field
- y component of the electric field
- A vector of six components of the electric field from dipole radiation
"""
function calEfield(ϕ, Tp, Ts, sinθ₁, cosθ₁; dvec = [1,1,1])

    pvec = [cosθ₁*cos(ϕ),cosθ₁*sin(ϕ),-sinθ₁].*Tp.*dvec
    svec = [-sin(ϕ),cos(ϕ),0.0].*Ts.*dvec
    hx = pvec.*cos(ϕ)-svec.*sin(ϕ)
    hy = pvec.*sin(ϕ)+svec.*cos(ϕ)
    return sum(hx), sum(hy), hcat(hx,hy)
end
