## Define types used in MicroscopePSFs
"""
`PSF` defines an abstract type
"""
abstract type PSF end


"""
    pdf(p::PSF,x_pixel::Tuple,x_emitter::Tuple)

abstract function: return the psf at pixel `x_pixel` with emitter located at `x_emitter`. 

# Arguments
- `p::PSF`              : psf structure
- `x_pixel::Tuple`      : location of pixel 
- `x_emitter::Tuple`    : location of emitter 


"""
function pdf(p::PSF,x_pixel::Tuple,x_emitter::Tuple)
end 

"""
    pdfₐ(p::PSF,x_pixel::Tuple,x_emitter::Tuple)

abstract function: return the complex amplitude at pixel `x_pixel` with emitter located at `x_emitter`. 

# Arguments
- `p::PSF`              : psf structure
- `x_pixel::Tuple`      : location of pixel 
- `x_emitter::Tuple`    : location of emitter 

"""
function pdfₐ(p::PSF,x_pixel::Tuple,x_emitter::Tuple)
end 

"""
    pdf(p::PSF,roi::Array,x_emitter::Tuple)

return the psf at pixel locations defined by `roi` with emitter located at `x_emitter`. 

# Arguments
- `p::PSF`              : psf structure
- `roi::Array`          : array of tuples that define the pixel locations 
- `x_emitter::Tuple`    : location of emitter 

"""
function pdf(p::PSF,roi::Array,x_emitter::Tuple)
    return pdf.(p,roi,Ref(x_emitter))
end

"""
    pdfₐ(p::PSF,roi::Array,x_emitter::Tuple)

return the complex amplitude at pixel locations defined by `roi` with emitter located at `x_emitter`. 

# Arguments
- `p::PSF`              : psf structure
- `roi::Array`          : array of tuples that define the pixel locations 
- `x_emitter::Tuple`    : location of emitter 

"""

function pdfₐ(p::PSF,roi::Array,x_emitter::Tuple)
    return pdfₐ.(p,roi,Ref(x_emitter))
end

"""
    pdf!(im::Array,p::PSF,roi::Array,x_emitter::Tuple)

update `im` to the psf at pixel locations defined by `roi` with emitter located at `x_emitter`. 

# Arguments
- `im::Array`           : array of numerical values with the same size as roi
- `p::PSF`              : psf structure
- `roi::Array`          : array of tuples that define the pixel locations 
- `x_emitter::Tuple`    : location of emitter 

"""

function pdf!(im::Array,p::PSF,roi::Array,x_emitter::Tuple)
    im.= pdf.(p,roi,Ref(x_emitter))
    return nothing
end

"""
    pdfₐ!(im::Array,p::PSF,roi::Array,x_emitter::Tuple)

update `im` to the complex amplitude at pixel locations defined by `roi` with emitter located at `x_emitter`. 

# Arguments
- `im::Array`           : array of numerical values with the same size as roi
- `p::PSF`              : psf structure
- `roi::Array`          : array of tuples that define the pixel locations 
- `x_emitter::Tuple`    : location of emitter 

"""

function pdfₐ!(im::Array,p::PSF,roi::Array,x_emitter::Tuple)
    im.= pdf.(p,roi,Ref(x_emitter))
    return nothing
end

"""
    pdf(p::PSF,roi::Array,x_emitter::Array)

return the overlap of multiple psfs at pixel locations defined by `roi` with emitter positions defined by `x_emitter`. 

# Arguments
- `p::PSF`              : psf structure
- `roi::Array`          : array of tuples that define the pixel locations 
- `x_emitter::Array`    : array of tuples that defines the emitter locations

"""

function pdf(p::PSF,roi::Array,x_emitter::Array)
    out=pdf.(p,roi,Ref(x_emitter[1]))
    for ii=2:length(x_emitter)
        out=out.+pdf.(p,roi,Ref(x_emitter[ii]))
    end
    return out
end

#needed for broadcasting
Base.broadcastable(x::PSF) = Ref(x)



