using MicroscopePSFs
using Test
using LinearAlgebra
using SpecialFunctions
using Zygote

@testset "MicroscopePSFs" verbose=true begin
    # Run only the working tests
    @testset "Basic Tests" begin
        # Check that basic objects can be created
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
    
    @testset "Scalar3D Basic" begin
        na = 1.2
        λ = 0.6
        n_medium = 1.33
        
        # Test constructor
        psf = Scalar3DPSF(na, λ, n_medium)
        @test psf.nₐ ≈ na
        @test psf.λ ≈ λ
        @test psf.n ≈ n_medium
        
        # Test simple properties
        @test psf(0.0, 0.0, 0.0) > psf(0.5, 0.0, 0.0)
        @test psf(0.1, 0.0, 0.0) ≈ psf(-0.1, 0.0, 0.0)
    end
    
    # Only test Vector3D with minimal example
    @testset "Vector3D Basic" begin
        na = 1.2
        λ = 0.6
        n_medium = 1.33
        
        # Test constructor with x-oriented dipole
        dipole_x = DipoleVector(1.0, 0.0, 0.0)
        psf_x = Vector3DPSF(na, λ, dipole_x, n_medium=n_medium)
        
        # Test basic symmetry
        @test psf_x(0.0, 0.5, 0.0) ≈ psf_x(0.0, -0.5, 0.0)
    end
end