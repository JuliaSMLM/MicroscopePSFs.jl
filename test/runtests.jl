using MicroscopePSFs
using Test

@testset "MicroscopePSFs.jl" begin
    PSF=MicroscopePSFs

    na=1.2
    n=1.3
    λ=.6 
    pixelsize=.1
    psf_airy=PSF.Airy2D(na,λ,pixelsize)
    psf_scalar=PSF.Scalar3D(na,λ,n,pixelsize)

    sz=4
    roi=[(i,j,0) for i=-sz/2:(sz/2-1), 
        j=-sz/2:(sz/2-1)] 

    im_airy=PSF.pdf(psf_airy,roi,(0.0,0.0,0.0))
    im_scalar=PSF.pdf(psf_scalar,roi,(0.0,0.0,0.0))
      
    @test minimum(isapprox(im_airy, im_scalar,atol=1e-4))

    ip=PSF.InterpolatedPSF(psf_scalar,(sz*2,sz*2,.2);subsampling=2)

    im_scalar=PSF.pdf(psf_scalar,roi,(sz/2+1,sz/2,.2))
    im_itp=PSF.pdf(ip,roi,(sz/2+1,sz/2,.2))

    @test minimum(isapprox(im_scalar, im_itp,atol=1e-4))



end
