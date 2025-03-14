using MicroscopePSFs
using Test
using LinearAlgebra
using SpecialFunctions
using HDF5

@testset "MicroscopePSFs" verbose=true begin
    
    # Basic PSFs Test: Gaussian2D and Airy2D
    @testset "Basic PSFs" begin
        include("basic_psfs_tests.jl")
    end
    
    # Scalar3D PSF tests
    @testset "Scalar3D PSF" begin
        include("scalar3d_tests.jl")
    end
    
    # Vector3D PSF tests
    @testset "Vector3D PSF" begin
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
    
    # I/O functionality tests
    @testset "I/O Functions" begin
        include("io_tests.jl")
    end
end