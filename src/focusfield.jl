

mutable struct Focusfield{PF<:PupilFunction,T<:AbstractFloat} <: PSF
    pupilfunctionx::PF
    pupilfunctiony::PF
    pupilfunctionz::PF
    pixelsize::T
    ksize::Int
    excitationfield::Vector
    normf::T
    electricfield::Char
    phasemask
end


function Focusfield(nₐ, λ, n::Vector, pixelsize;
    normf=1.0,
    ksize=256,
    electricfield='x',
    excitationfield=[1.0, 0],
    fθ=(x) -> 0.0,
    fϕ=(x) -> 0.0,
    zstage=0.0,
    sigmax_pixel= 0.1 # in pixel
)

    pupilx = zeros(ksize, ksize, 2)
    pupily = zeros(ksize, ksize, 2)
    pupilz = zeros(ksize, ksize, 2)
    phasemask = zeros(ksize, ksize)
    kpixelsize = 2 * nₐ / λ / ksize
    k = n[1] / λ
    k0 = (ksize + 1) / 2

    sigmax = pixelsize *sigmax_pixel
    sigmak = 1 / (2 * pi * sigmax)
    println("sigmak=$sigmak")
    println("kpixelsize=$kpixelsize")


    for ii in 1:ksize, jj in 1:ksize # jj is index along y  
        kx = kpixelsize * (ii - k0)
        ky = kpixelsize * (jj - k0)
        kr2 = kx^2 + ky^2

        Tp, Ts, sinθ₁, cosθ₁, _, cosθ₃ = calFresnel(kr2, λ, n)
        immphase = exp(-2*pi*(n[1]/λ*cosθ₁*zstage)*im)

        if kr2 < (nₐ / λ)^2
            g = exp(-kr2/(2*sigmak^2))
            ϕ = atan(ky, kx)
            θ = asin(sqrt(kr2) / k)
            nθ = [cos(ϕ), sin(ϕ)]
            nϕ = [-sin(ϕ), cos(ϕ)]
            Eθ = dot(nθ, excitationfield.*g)
            Eϕ = dot(nϕ, excitationfield.*g)
            #Eθ = 0.0
            #Eϕ = 1.0
            # transmit through sample medium
            Eθt = Eθ * Tp * exp(im * fθ(θ) + im * fϕ(ϕ))*immphase
            Eϕt = Eϕ * Ts * exp(im * fθ(θ) + im * fϕ(ϕ))*immphase
            # decompose into Ex, Ey, Ez
            Ez = Eθt * sinθ₁ * n[1] / n[3]
            Ex = Eθt * cosθ₃ * cos(ϕ) - Eϕt * sin(ϕ)
            Ey = Eθt * cosθ₃ * sin(ϕ) + Eϕt * cos(ϕ)

            phasemask[jj, ii] = fθ(θ) + fϕ(ϕ)
            pupilx[jj, ii, 1] = abs(Ex) * kpixelsize
            pupilx[jj, ii, 2] = angle(Ex)
            pupily[jj, ii, 1] = abs(Ey) * kpixelsize
            pupily[jj, ii, 2] = angle(Ey)
            pupilz[jj, ii, 1] = abs(Ez) * kpixelsize
            pupilz[jj, ii, 2] = angle(Ez)
        end
    end

    px = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupilx)
    py = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupily)
    pz = PupilFunction(nₐ, λ, n[3], pixelsize, kpixelsize, pupilz)

    return Focusfield{PupilFunction,typeof(nₐ)}(px, py, pz, pixelsize, ksize, excitationfield, normf, electricfield,phasemask)

end


function pdfₐ(p::Focusfield, pixel::Tuple, x_emitter::Tuple)
    if p.electricfield == 'x'
        out = pdfₐ(p.pupilfunctionx, pixel, x_emitter) * p.pixelsize / p.normf
    end
    if p.electricfield == 'y'
        out = pdfₐ(p.pupilfunctiony, pixel, x_emitter) * p.pixelsize / p.normf
    end
    if p.electricfield == 'z'
        out = pdfₐ(p.pupilfunctionz, pixel, x_emitter) * p.pixelsize / p.normf
    end
    return out
end

function pdfₐ(p::Focusfield, roi::Array, x_emitter::Tuple)
    out = zeros(Complex, size(roi))
    Threads.@threads for ii = 1:length(roi)
        out[ii] = pdfₐ(p, roi[ii], x_emitter)
    end
    return out
end

function pdf(p::Focusfield, pixel::Tuple, x_emitter::Tuple)
    return abs2(pdfₐ(p, pixel, x_emitter))
end

function pdf(p::Focusfield, roi::Array, x_emitter::Tuple)
    out = zeros(size(roi))
    Threads.@threads for ii = 1:length(roi)
        out[ii] = pdf(p, roi[ii], x_emitter)
    end
    return out
end