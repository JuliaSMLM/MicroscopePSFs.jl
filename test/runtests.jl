using MicroscopePSFs
using Test
using LinearAlgebra
using SpecialFunctions
using HDF5

@testset "MicroscopePSFs" verbose=true begin
    
    # Basic PSFs Test: GaussianPSF and AiryPSF
    @testset "Basic PSFs" begin
        include("basic_psfs_tests.jl")
    end
    
    # ScalarPSF tests
    @testset "ScalarPSF" begin
        include("scalar3d_tests.jl")
    end
    
    # VectorPSF tests
    @testset "VectorPSF" begin
        include("vector3d_tests.jl")
    end
    
    # Spline PSF tests
    @testset "Spline PSF" begin
        include("spline_psf_tests.jl")
    end
    
    # Aberration and Zernike tests
    @testset "Aberrations" begin
        include("aberration_tests.jl")
    end

     # Multi-emitter tests
     @testset "Multi Emitter" begin
        include("multi_emitter_tests.jl")
    end

    # Support region tests
    @testset "Support Region" begin
        include("support_region_tests.jl")
    end
    
    # I/O functionality tests
    @testset "I/O Functions" begin
        include("io_tests.jl")
    end

   
end