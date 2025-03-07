# test/scalar3d_tests.jl
using Test
using MicroscopePSFs

@testset "Scalar3D PSF" begin
    # Setup common parameters
    na = 1.2
    λ = 0.6
    n_medium = 1.33
    
    @testset "Constructor" begin
        # Test basic constructor
        psf = Scalar3DPSF(na, λ, n_medium)
        @test psf.nₐ ≈ na
        @test psf.λ ≈ λ
        @test psf.n ≈ n_medium
        
        # Test with Zernike coefficients
        zernike = ZernikeCoefficients([0.1, 0.2, 0.0, 0.3])
        psf_zernike = Scalar3DPSF(na, λ, n_medium; coeffs=zernike)
        @test psf_zernike.zernike_coeffs.coeffs[2] ≈ 0.2
    end
    
    @testset "Physical Correctness" begin
        psf = Scalar3DPSF(na, λ, n_medium)
        
        # Test central maximum
        @test psf(0.0, 0.0, 0.0) > psf(0.1, 0.0, 0.0)
        @test psf(0.0, 0.0, 0.0) > psf(0.0, 0.1, 0.0)
        
        # Test symmetry
        @test psf(0.1, 0.0, 0.0) ≈ psf(-0.1, 0.0, 0.0)
        @test psf(0.0, 0.1, 0.0) ≈ psf(0.0, -0.1, 0.0)
        @test psf(0.1, 0.1, 0.0) ≈ psf(-0.1, -0.1, 0.0)
        
        # Test z-dependence (should have lower intensity away from focus)
        @test psf(0.0, 0.0, 0.0) > psf(0.0, 0.0, 1.0)
        
        # Test that far from PSF center, intensity approaches zero
        @test psf(10.0, 0.0, 0.0) < 1e-3
    end
    
    @testset "Integration" begin
        psf = Scalar3DPSF(na, λ, n_medium)
        
        # Create camera and emitter
        camera = IdealCamera(0:0.1:2.0, 0:0.1:2.0)
        emitter = Emitter3D(1.0, 1.0, 0.0, 1000.0)
        
        # Integrate PSF over pixels
        img = integrate_pixels(psf, camera, emitter)
        
        # Check image properties
        @test size(img) == (length(camera.y_pixels), length(camera.x_pixels))
        
        # Total photons should be close to emitter photons (energy conservation)
        @test sum(img) > 0.9 * emitter.photons
        
        # Peak should be near emitter position
        max_idx = findmax(img)[2]
        center_x = findfirst(x -> x ≈ 1.0, camera.x_pixels)
        center_y = findfirst(y -> y ≈ 1.0, camera.y_pixels)
        @test abs(max_idx[1] - center_y) ≤ 1
        @test abs(max_idx[2] - center_x) ≤ 1
    end
    
    @testset "Aberrations" begin
        # Create base PSF
        psf = Scalar3DPSF(na, λ, n_medium)
        
        # Create aberrated PSF with astigmatism
        zernike = ZernikeCoefficients(zeros(5))
        zernike.coeffs[5] = 0.2  # Astigmatism
        psf_aberrated = Scalar3DPSF(na, λ, n_medium; coeffs=zernike)
        
        # Test that on-axis points should still match
        @test psf(0.0, 0.0, 0.0) ≈ psf_aberrated(0.0, 0.0, 0.0) rtol=1e-3
        
        # Astigmatism breaks circular symmetry
        # Test at points on the x and y axes (should be different)
        @test abs(psf_aberrated(0.5, 0.0, 0.0) - psf_aberrated(0.0, 0.5, 0.0)) > 1e-3
    end
end