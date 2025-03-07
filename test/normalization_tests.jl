# test/normalization_tests.jl
using Test
using MicroscopePSFs

@testset "Normalization" begin
    @testset "Energy Conservation - Scalar3D" begin
        # Create PSF
        na = 1.2
        λ = 0.6
        n_medium = 1.33
        psf = Scalar3DPSF(na, λ, n_medium)
        
        # Create camera with reasonable field of view
        camera = IdealCamera(-1.0:0.1:1.0, -1.0:0.1:1.0)
        
        # Test energy conservation
        photons = 1000.0
        emitter = Emitter3D(0.0, 0.0, 0.0, photons)
        
        # Generate image
        img = integrate_pixels(psf, camera, emitter)
        
        # Total energy should be within reasonable range of emitter photons
        # Some energy can be lost due to finite image size
        @test sum(img) > 0.8 * photons
        @test sum(img) < 1.1 * photons
    end
end