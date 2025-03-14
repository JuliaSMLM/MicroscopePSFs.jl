# test/aberration_tests.jl

@testset "Zernike Coefficients" begin
    # Test ZernikeCoefficients constructor
    mag = [0.1, 0.2, 0.3, 0.4, 0.5]
    phase = zeros(5)  # Create matching phase vector
    zc = ZernikeCoefficients(mag, phase)
    
    # Test index access with correct field names
    @test zc.mag[1] ≈ 0.1
    @test zc.mag[5] ≈ 0.5
    
    # Test length
    @test length(zc) == 5
    
    # Test zero initialization
    zc_zeros = ZernikeCoefficients(5)
    # Note: default constructor sets mag[1] = 1.0, not 0
    @test zc_zeros.mag[1] ≈ 1.0
    @test all(zc_zeros.mag[2:end] .== 0.0)
    @test all(zc_zeros.phase .== 0.0)
end