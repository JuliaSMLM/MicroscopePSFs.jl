# test/psf_base_tests.jl
using Test
using MicroscopePSFs
using LinearAlgebra

@testset "PSF Base Types" begin
    # Check that basic objects can be created
    @testset "Basic Object Creation" begin
        @test isa(ZernikeCoefficients(5), ZernikeCoefficients)
        @test isa(Gaussian2D(0.15), Gaussian2D)
        @test isa(DipoleVector(1.0, 0.0, 0.0), DipoleVector)
        
        # Test Emitter types
        emitter2d = Emitter2D(1.5, 2.0, 1000.0)
        @test emitter2d.x ≈ 1.5
        @test emitter2d.y ≈ 2.0
        @test emitter2d.photons ≈ 1000.0
        
        emitter3d = Emitter3D(1.5, 2.0, 3.0, 1000.0)
        @test emitter3d.x ≈ 1.5
        @test emitter3d.y ≈ 2.0
        @test emitter3d.z ≈ 3.0
        @test emitter3d.photons ≈ 1000.0
        
        # Test basic PSF evaluation
        psf_gauss = Gaussian2D(0.15)
        @test psf_gauss(0.0, 0.0) > psf_gauss(0.1, 0.1)
    end

    # Test basic PSF types
    @testset "2D PSF Constructors" begin
        # Test Airy PSF constructor
        na = 1.2
        λ = 0.6
        pixelsize = 0.1
        
        psf_airy = Airy2D(na, λ, pixelsize)
        @test psf_airy.na ≈ na
        @test psf_airy.λ ≈ λ
        @test psf_airy.pixelsize ≈ pixelsize
        
        # Test Gaussian PSF constructor
        sigma = 0.15
        psf_gauss = Gaussian2D(sigma)
        @test psf_gauss.sigma ≈ sigma
        
        # Test physical reasonableness
        xy_points = [(0.0, 0.0), (0.1, 0.0), (0.0, 0.1), (0.1, 0.1)]
        
        # PSF should be maximum at center
        airy_vals = [psf_airy(x, y) for (x, y) in xy_points]
        @test argmax(airy_vals) == 1
        
        # PSF should decrease with radius
        @test airy_vals[1] > airy_vals[2]
        @test airy_vals[1] > airy_vals[3]
        @test airy_vals[1] > airy_vals[4]
    end
    
    @testset "Emitter Types" begin
        # Test emitter constructors
        emitter2d = Emitter2D(1.5, 2.0, 1000.0)
        @test emitter2d.x ≈ 1.5
        @test emitter2d.y ≈ 2.0
        @test emitter2d.photons ≈ 1000.0
        
        emitter3d = Emitter3D(1.5, 2.0, 3.0, 1000.0)
        @test emitter3d.x ≈ 1.5
        @test emitter3d.y ≈ 2.0
        @test emitter3d.z ≈ 3.0
        @test emitter3d.photons ≈ 1000.0
        
        # Test dipole emitter
        dipole_emitter = DipoleEmitter3D(1.5, 2.0, 3.0, 1000.0, 1.0, 0.0, 0.0)
        @test dipole_emitter.x ≈ 1.5
        @test dipole_emitter.y ≈ 2.0
        @test dipole_emitter.z ≈ 3.0
        @test dipole_emitter.photons ≈ 1000.0
        
        # Test DipoleVector
        dipole = DipoleVector(1.0, 2.0, 3.0)
        @test dipole.px^2 + dipole.py^2 + dipole.pz^2 ≈ 1.0  # Should be normalized
        
        # Test that unnormalized vector gets normalized
        v = [3.0, 0.0, 4.0]  # 3-4-5 triangle
        v_norm = v / norm(v)
        dipole = DipoleVector(v[1], v[2], v[3])
        @test dipole.px ≈ v_norm[1]
        @test dipole.py ≈ v_norm[2]
        @test dipole.pz ≈ v_norm[3]
    end
    
    @testset "Camera Types" begin
        # Test camera constructor
        x_pixels = 0:0.1:1.0
        y_pixels = 0:0.1:1.0
        camera = IdealCamera(x_pixels, y_pixels)
        
        @test camera.x_pixels == x_pixels
        @test camera.y_pixels == y_pixels
        @test size(camera) == (length(y_pixels), length(x_pixels))
        
        # Test camera bounds
        @test minimum(camera.x_pixels) ≈ 0.0
        @test maximum(camera.x_pixels) ≈ 1.0
        @test minimum(camera.y_pixels) ≈ 0.0
        @test maximum(camera.y_pixels) ≈ 1.0
    end
end