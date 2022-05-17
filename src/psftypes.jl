## Define types used in MicroscopePSFs
"""
    `PSF` defines an abstract type
"""
abstract type PSF end


"""
    pdf(p::PSF,x_pixel::Tuple,x_emitter::Tuple)

return the psf at pixel `x_pixel` with emitter located at `x_emitter`. 

#Arguments
- `p::PSF`              : psf structure
- `x_pixel::Tuple`      : location of pixel 
- `x_emitter::Tuple`    : location of emitter 


"""
function pdf(p::PSF,x_pixel::Tuple,x_emitter::Tuple)
end 

"""
    pdfₐ(p::PSF,x_pixel::Tuple,x_emitter::Tuple)

return the complex amplitude at pixel `x_pixel` with emitter located at `x_emitter`. 

#Arguments
- `p::PSF`              : psf structure
- `x_pixel::Tuple`      : location of pixel 
- `x_emitter::Tuple`    : location of emitter 

"""
function pdfₐ(p::PSF,x_pixel::Tuple,x_emitter::Tuple)
end 


function pdf(p::PSF,roi::Array,x_emitter::Tuple)
    return pdf.(p,roi,Ref(x_emitter))
end

function pdfₐ(p::PSF,roi::Array,x_emitter::Tuple)
    return pdfₐ.(p,roi,Ref(x_emitter))
end

function pdf!(im::Array,p::PSF,roi::Array,x_emitter::Tuple)
    im.= pdf.(p,roi,Ref(x_emitter))
    return nothing
end

function pdfₐ!(im::Array,p::PSF,roi::Array,x_emitter::Tuple)
    im.= pdf.(p,roi,Ref(x_emitter))
    return nothing
end

function pdf(p::PSF,roi::Array,x_emitter::Array)
    out=pdf.(p,roi,Ref(x_emitter[1]))
    for ii=2:length(x_emitter)
        out=out.+pdf.(p,roi,Ref(x_emitter[ii]))
    end
    return out
end

#needed for broadcasting
Base.broadcastable(x::PSF) = Ref(x)



