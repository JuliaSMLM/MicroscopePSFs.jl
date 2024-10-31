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
mutable struct Dipole3D_fast{PF<:PupilFunction1d,T<:AbstractFloat,I<:Int} <: PSF
    pupilfunctionx::PF
    pupilfunctiony::PF 
    pixelsize::T
    Σ::T
    dipole_ang::Vector
    ksize::Int
    electricfield::Char
    normf::T
    excitationfield::Union{Vector,Float64}
end


function Dipole3D_fast(nₐ, λ, n::Vector, pixelsize, dipole_ang::Vector; 
    normf=1.0, 
    zstage=0.0, 
    excitationfield=1.0,
    electricfield='x', 
    Σ=0.0, 
    ksize=256, 
    mvtype="bead",
    δ::Union{Nothing,Float64}=nothing)

    pupilx = zeros(1,ksize)
    pupily = zeros(1,ksize)
    #apod = zeros(ksize,ksize,2)
    kpixelsize = nₐ / λ / ksize
    
    α = dipole_ang[1]
    β = dipole_ang[2]
    dvec = [sin(α)*cos(β),sin(α)*sin(β),cos(α)]

    gx,gy = getfield_fun(n,λ,zstage,excitationfield,δ; dvec=dvec)


    if mvtype == "bead"
        px = PupilFunction1d(nₐ, λ, n[1], pixelsize, kpixelsize, pupilx,gx)
        py = PupilFunction1d(nₐ, λ, n[1], pixelsize, kpixelsize, pupily,gy)
    elseif mvtype == "stage"
        px = PupilFunction1d(nₐ, λ, n[3], pixelsize, kpixelsize, pupilx,gx)
        py = PupilFunction1d(nₐ, λ, n[3], pixelsize, kpixelsize, pupily,gy)
    else
        error("mvtype not recognized")
    end
    #println(typeof(px))
    #normalize!(px)
    #normalize!(py)
    return Dipole3D_fast{PupilFunction1d,typeof(nₐ),Int}(px,py, pixelsize, Σ,dipole_ang,ksize,electricfield,normf,excitationfield)
end

function pdfₐ(p::Dipole3D_fast, pixel::Tuple,x_emitter::Tuple)
    if p.electricfield == 'x'
        out = pdfₐ(p.pupilfunctionx, pixel, x_emitter) * p.pixelsize / p.normf
    end
    if p.electricfield == 'y'
        out = pdfₐ(p.pupilfunctiony, pixel, x_emitter) * p.pixelsize / p.normf
    end
    return out
end    

function pdfₐ(p::Dipole3D_fast, roi::Array,x_emitter::Tuple)
    out=zeros(Complex,size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdfₐ(p, roi[ii], x_emitter)
    end
    return out
end    

function pdf(p::Dipole3D_fast, pixel::Tuple,x_emitter::Tuple)
    out = pdfₐ(p, pixel, x_emitter)
    return real(out*conj(out))
end    

function pdf(p::Dipole3D_fast,roi::Array,x_emitter::Tuple)
    out=zeros(size(roi))
    Threads.@threads for ii=1:length(roi)
        out[ii]=pdf(p,roi[ii],x_emitter)
    end
    return out
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
function calEfield_1d(kr,x,y,n,λ,zstage,Eex::Vector,δ::Nothing; dvec = [1,1,1])
    #dvec = dvec
    #Eex = Eex
    normdvec = norm(dvec)
    normEex = norm(Eex)
    E0 = dot(Eex,dvec)/normdvec/normEex/normdvec
    #dvec = dvec.*E0
    kr2=kr^2
    Tp, Ts, sinθ₁, cosθ₁,  _, cosθ₃ = calFresnel(kr2,λ,n)
    immphase = exp(-2*pi*(n[3]/λ*cosθ₃*zstage)*im)
    apod = sqrt(cosθ₃)/cosθ₁
    cos2ϕᵣ,sin2ϕᵣ,cosϕᵣ, sinϕᵣ, r = car2pol(x,y)
    J0 = bJ(0,kr*r)
    J1 = bJ(1,kr*r)
    J2 = bJ(2,kr*r)

    #hx = [cosθ₁*pi*(J0-J2*cos2ϕᵣ),-cosθ₁*pi*J2*sin2ϕᵣ,-sinθ₁*2*pi*im*J1*cosϕᵣ].*Tp.*dvec .-
    #    [-pi*(J0+J2*cos2ϕᵣ),-pi*J2*sin2ϕᵣ,0.0].*Ts.*dvec

    #hy = [-cosθ₁*pi*J2*sin2ϕᵣ,cosθ₁*pi*(J0+J2*cos2ϕᵣ),-sinθ₁*2*pi*im*J1*sinϕᵣ].*Tp.*dvec .+
    #    [pi*J2*sin2ϕᵣ,pi*(J0-J2*cos2ϕᵣ),0.0].*Ts.*dvec

    hx = cosθ₁*pi*(J0-J2*cos2ϕᵣ)*Tp*dvec[1]-cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[2]-sinθ₁*2*pi*im*J1*cosϕᵣ*Tp*dvec[3]+
        pi*(J0+J2*cos2ϕᵣ)*Ts*dvec[1]+pi*J2*sin2ϕᵣ*Ts*dvec[2]

    hy = -cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[1]+cosθ₁*pi*(J0+J2*cos2ϕᵣ)*Tp*dvec[2]-sinθ₁*2*pi*im*J1*sinϕᵣ*Tp*dvec[3]+
        pi*J2*sin2ϕᵣ*Ts*dvec[1]+pi*(J0-J2*cos2ϕᵣ)*Ts*dvec[2]
        
    hx = hx*immphase*apod*E0
    hy = hy*immphase*apod*E0

    return hx, hy
end

function getfield_fun(n,λ,zstage,Eex::Union{Vector,Float64},δ::Union{Nothing,Float64}; dvec = [1,1,1])
    fx = (kr,x,y) -> calEfield_1d(kr,x,y,n,λ,zstage,Eex, δ; dvec=dvec)[1]
    fy = (kr,x,y) -> calEfield_1d(kr,x,y,n,λ,zstage,Eex, δ; dvec=dvec)[2]
    return fx,fy
end


function calEfield_1d(kr,x,y,n,λ,zstage,Eex::Float64,δ::Nothing; dvec = [1,1,1])
    #dvec = dvec./norm(dvec)
    #dvec = dvec.*Eex
    normdvec = norm(dvec)
    kr2=kr^2
    Tp, Ts, sinθ₁, cosθ₁,  _, cosθ₃ = calFresnel(kr2,λ,n)
    immphase = exp(-2*pi*(n[3]/λ*cosθ₃*zstage)*im)
    apod = sqrt(cosθ₃)/cosθ₁
    cos2ϕᵣ,sin2ϕᵣ,cosϕᵣ, sinϕᵣ, r = car2pol(x,y)
    J0 = bJ(0,kr*r)
    J1 = bJ(1,kr*r)
    J2 = bJ(2,kr*r)

    #hx = [cosθ₁*pi*(J0-J2*cos2ϕᵣ),-cosθ₁*pi*J2*sin2ϕᵣ,-sinθ₁*2*pi*im*J1*cosϕᵣ].*Tp.*dvec .-
    #    [-pi*(J0+J2*cos2ϕᵣ),-pi*J2*sin2ϕᵣ,0.0].*Ts.*dvec


    hx = cosθ₁*pi*(J0-J2*cos2ϕᵣ)*Tp*dvec[1]-cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[2]-sinθ₁*2*pi*im*J1*cosϕᵣ*Tp*dvec[3]+
        pi*(J0+J2*cos2ϕᵣ)*Ts*dvec[1]+pi*J2*sin2ϕᵣ*Ts*dvec[2]
    

    #hy = [-cosθ₁*pi*J2*sin2ϕᵣ,cosθ₁*pi*(J0+J2*cos2ϕᵣ),-sinθ₁*2*pi*im*J1*sinϕᵣ].*Tp.*dvec .+
    #    [pi*J2*sin2ϕᵣ,pi*(J0-J2*cos2ϕᵣ),0.0].*Ts.*dvec

    
    hy = -cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[1]+cosθ₁*pi*(J0+J2*cos2ϕᵣ)*Tp*dvec[2]-sinθ₁*2*pi*im*J1*sinϕᵣ*Tp*dvec[3]+
        pi*J2*sin2ϕᵣ*Ts*dvec[1]+pi*(J0-J2*cos2ϕᵣ)*Ts*dvec[2]

    hx = hx*immphase*apod*Eex/normdvec
    hy = hy*immphase*apod*Eex/normdvec
    return hx, hy
end

function calEfield_1d(kr,x,y,n,λ,zstage,Eex::Float64,δ::Float64; dvec = [1,1,1])
    #dvec = dvec./norm(dvec)
    #dvec = dvec.*Eex
    normdvec = norm(dvec)
    kr2=kr^2
    Tp, Ts, sinθ₁, cosθ₁,  _, cosθ₃ = calFresnel(kr2,λ,n)
    immphase = exp(-2*pi*(n[3]/λ*cosθ₃*zstage)*im)
    apod = sqrt(cosθ₃)/cosθ₁
    cos2ϕᵣ,sin2ϕᵣ,cosϕᵣ, sinϕᵣ, r = car2pol(x,y)
    J0 = bJ(0,kr*r)
    J1 = bJ(1,kr*r)
    J2 = bJ(2,kr*r)

    hx = cosθ₁*pi*(J0-J2*cos2ϕᵣ)*Tp*dvec[1]-cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[2]-sinθ₁*2*pi*im*J1*cosϕᵣ*Tp*dvec[3]+
        pi*(J0+J2*cos2ϕᵣ)*Ts*dvec[1]+pi*J2*sin2ϕᵣ*Ts*dvec[2]

    hy = -cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[1]+cosθ₁*pi*(J0+J2*cos2ϕᵣ)*Tp*dvec[2]-sinθ₁*2*pi*im*J1*sinϕᵣ*Tp*dvec[3]+
        pi*J2*sin2ϕᵣ*Ts*dvec[1]+pi*(J0-J2*cos2ϕᵣ)*Ts*dvec[2]

    #h_rad = [cosθ₁*cosϕᵣ*2*pi*im*J1,cosθ₁*sinϕᵣ*2*pi*im*J1,-sinθ₁*2*pi*J0].*Tp.*dvec
    #h_azi = [-sinϕᵣ*2*pi*im*J1,cosϕᵣ*2*pi*im*J1,0.0].*Ts.*dvec

    h_rad = cosθ₁*cosϕᵣ*2*pi*im*J1*Tp*dvec[1]+cosθ₁*sinϕᵣ*2*pi*im*J1*Tp*dvec[2]-sinθ₁*2*pi*J0*Tp*dvec[3]
    h_azi = -sinϕᵣ*2*pi*im*J1*Ts*dvec[1]+cosϕᵣ*2*pi*im*J1*Ts*dvec[2]

    h_qx = h_rad*cos(δ/2)*im-hx*sin(δ/2)
    h_qy = -h_azi*cos(δ/2)*im-hy*sin(δ/2)

    h_qx = h_qx*immphase*apod*Eex/normdvec
    h_qy = h_qy*immphase*apod*Eex/normdvec

    return h_qx,h_qy
end

function calEfield_1d(kr,x,y,n,λ,zstage,Eex::Vector,δ::Float64; dvec = [1,1,1])
    #dvec = dvec./norm(dvec)
    #Eex = Eex./norm(Eex)
    normdvec = norm(dvec)
    normEex = norm(Eex)
    E0 = dot(Eex,dvec)/normdvec/normEex/normdvec
    #dvec = dvec.*E0
    kr2=kr^2
    Tp, Ts, sinθ₁, cosθ₁,  _, cosθ₃ = calFresnel(kr2,λ,n)
    immphase = exp(-2*pi*(n[3]/λ*cosθ₃*zstage)*im)
    apod = sqrt(cosθ₃)/cosθ₁
    cos2ϕᵣ,sin2ϕᵣ,cosϕᵣ, sinϕᵣ, r = car2pol(x,y)
    J0 = bJ(0,kr*r)
    J1 = bJ(1,kr*r)
    J2 = bJ(2,kr*r)

    hx = cosθ₁*pi*(J0-J2*cos2ϕᵣ)*Tp*dvec[1]-cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[2]-sinθ₁*2*pi*im*J1*cosϕᵣ*Tp*dvec[3]+
        pi*(J0+J2*cos2ϕᵣ)*Ts*dvec[1]+pi*J2*sin2ϕᵣ*Ts*dvec[2]

    hy = -cosθ₁*pi*J2*sin2ϕᵣ*Tp*dvec[1]+cosθ₁*pi*(J0+J2*cos2ϕᵣ)*Tp*dvec[2]-sinθ₁*2*pi*im*J1*sinϕᵣ*Tp*dvec[3]+
        pi*J2*sin2ϕᵣ*Ts*dvec[1]+pi*(J0-J2*cos2ϕᵣ)*Ts*dvec[2]

    h_rad = cosθ₁*cosϕᵣ*2*pi*im*J1*Tp*dvec[1]+cosθ₁*sinϕᵣ*2*pi*im*J1*Tp*dvec[2]-sinθ₁*2*pi*J0*Tp*dvec[3]
    h_azi = -sinϕᵣ*2*pi*im*J1*Ts*dvec[1]+cosϕᵣ*2*pi*im*J1*Ts*dvec[2]

    h_qx = h_rad*cos(δ/2)*im-hx*sin(δ/2)
    h_qy = -h_azi*cos(δ/2)*im-hy*sin(δ/2)

    h_qx = h_qx*immphase*apod*E0
    h_qy = h_qy*immphase*apod*E0

    return h_qx,h_qy
end



function car2pol(x,y)
    r = sqrt(x^2 + y^2)
    if r == 0
        return 0.0,0.0,1.0,0.0,0.0
    end
    sinϕᵣ = y/r
    cosϕᵣ = x/r
    cos2ϕᵣ = 2*cosϕᵣ^2-1
    sin2ϕᵣ = 2*sinϕᵣ*cosϕᵣ
    return cos2ϕᵣ,sin2ϕᵣ,cosϕᵣ, sinϕᵣ, r
end

function bJ(n, x)
    return SpecialFunctions.besselj(n,2*pi*x)
end