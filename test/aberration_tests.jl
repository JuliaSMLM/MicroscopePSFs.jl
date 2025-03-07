# test/aberration_tests.jl
using Test
using MicroscopePSFs

@testset "Aberrations" begin
    @testset "Zernike Coefficients" begin
        # Test ZernikeCoefficients constructor
        coeffs = [0.1, 0.2, 0.3, 0.4, 0.5]
        zc = ZernikeCoefficients(coeffs)
        
        # Test index access
        @test zc.coeffs[1] ≈ 0.1
        @test zc.coeffs[5] ≈ 0.5
        
        # Test length
        @test length(zc) == 5
        
        # Test zero initialization
        zc_zeros = ZernikeCoefficients(5)
        @test all(zc_zeros.coeffs .== 0.0)
    end
end