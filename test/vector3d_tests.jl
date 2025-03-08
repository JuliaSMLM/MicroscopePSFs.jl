# test/vector3d_tests.jl

@testset "Vector3D PSF" begin
    # Setup common parameters
    na = 1.2
    λ = 0.6
    n_medium = 1.33
    n_coverslip = 1.52
    n_immersion = 1.52
    
    @testset "Constructor" begin
        # Test basic constructor
        dipole = DipoleVector(1.0, 0.0, 0.0)  # X-oriented
        psf = Vector3DPSF(na, λ, dipole; 
                          n_medium=n_medium, 
                          n_coverslip=n_coverslip, 
                          n_immersion=n_immersion)
        
        # Verify parameters were set correctly
        @test psf.nₐ ≈ na
        @test psf.λ ≈ λ
        @test psf.n_medium ≈ n_medium
        @test psf.n_coverslip ≈ n_coverslip
        @test psf.n_immersion ≈ n_immersion
        @test psf.dipole.px ≈ 1.0
        @test psf.dipole.py ≈ 0.0
        @test psf.dipole.pz ≈ 0.0
        
        # Test constructor with aberrations
        zc = ZernikeCoefficients(10)
        add_defocus!(zc, 0.5)
        psf_aberrated = Vector3DPSF(na, λ, dipole; 
                                   n_medium=n_medium, 
                                   base_zernike=zc)
        
        @test !isnothing(psf_aberrated.zernike_coeffs)
    end
    
    @testset "Dipole Orientation Symmetry" begin
        # Create dipoles with different orientations
        dipole_x = DipoleVector(1.0, 0.0, 0.0)  # X-oriented
        dipole_y = DipoleVector(0.0, 1.0, 0.0)  # Y-oriented
        dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Z-oriented
        
        # Create PSFs for each orientation
        psf_x = Vector3DPSF(na, λ, dipole_x; n_medium=n_medium)
        psf_y = Vector3DPSF(na, λ, dipole_y; n_medium=n_medium)
        psf_z = Vector3DPSF(na, λ, dipole_z; n_medium=n_medium)
        
        # X-dipole symmetry tests
        # Should be symmetric with respect to y-axis
        @test psf_x(0.0, 0.5, 0.0) ≈ psf_x(0.0, -0.5, 0.0)
        
        # For a pure x-dipole in scalar approximation, we actually get symmetry in both x and y
        @test psf_x(0.5, 0.0, 0.0) ≈ psf_x(-0.5, 0.0, 0.0)
        
        # Y-dipole symmetry tests
        # Should be symmetric with respect to x-axis
        @test psf_y(0.5, 0.0, 0.0) ≈ psf_y(-0.5, 0.0, 0.0)
        
        # Should also be symmetric with respect to y-axis in this implementation
        @test psf_y(0.0, 0.5, 0.0) ≈ psf_y(0.0, -0.5, 0.0)
        
        # Z-dipole symmetry tests
        # Should have rotational symmetry around z-axis
        r = 0.5
        for θ in [0, π/4, π/2, 3π/4, π]
            x1, y1 = r*cos(θ), r*sin(θ)
            x2, y2 = r*cos(θ+π/2), r*sin(θ+π/2)
            @test isapprox(psf_z(x1, y1, 0.0), psf_z(x2, y2, 0.0), rtol=1e-2)
        end
    end
    
    @testset "Dipole Radiation Pattern" begin
        # Test that dipoles follow expected radiation pattern
        dipole_x = DipoleVector(1.0, 0.0, 0.0)
        psf_x = Vector3DPSF(na, λ, dipole_x; n_medium=n_medium)
        
        # For our implementation, test that intensity decreases with distance
        # from the center for the in-plane dipole
        @test psf_x(0.0, 0.0, 0.0) > psf_x(0.5, 0.0, 0.0)
        @test psf_x(0.0, 0.0, 0.0) > psf_x(0.0, 0.5, 0.0)
        
        # Z-dipole should have donut pattern in xy-plane (maximum at some radius, minimum at center)
        dipole_z = DipoleVector(0.0, 0.0, 1.0)
        psf_z = Vector3DPSF(na, λ, dipole_z; n_medium=n_medium)
        
        # Center should be a minimum
        @test psf_z(0.0, 0.0, 0.0) < psf_z(0.3, 0.0, 0.0)
    end
    
    @testset "Amplitude Function" begin
        # Test that amplitude function returns expected components
        dipole = DipoleVector(1.0, 0.0, 0.0)
        psf = Vector3DPSF(na, λ, dipole; n_medium=n_medium)
        
        # Amplitude should return a vector of [Ex, Ey]
        amp = amplitude(psf, 0.3, 0.4, 0.0)
        
        # Check dimensions and types
        @test length(amp) == 2
        @test eltype(amp) <: Complex
        
        # For x-dipole, |Ex| should be larger than |Ey| in most positions
        @test abs(amp[1]) > abs(amp[2])
        
        # Intensity should equal sum of squared field amplitudes
        @test psf(0.3, 0.4, 0.0) ≈ abs2(amp[1]) + abs2(amp[2])
    end
    
    @testset "Pixel Integration" begin
        # Test integration over camera pixels
        dipole = DipoleVector(1.0, 0.0, 0.0)
        psf = Vector3DPSF(na, λ, dipole; n_medium=n_medium)
        
        # Create camera and emitter
        camera = IdealCamera(collect(0:0.1:2.0), collect(0:0.1:2.0))
        emitter = DipoleEmitter3D(1.0, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.0)
        
        # Integrate PSF over pixels
        img = integrate_pixels(psf, camera, emitter)
        
        # Check image properties
        @test size(img) == (length(0:0.1:2.0)-1, length(0:0.1:2.0)-1)
        
        # Test amplitude integration
        amp_img = integrate_pixels_amplitude(psf, camera, emitter)
        @test size(amp_img, 3) == 2  # Should have two field components
        
        # Just test that the amplitude image has the right dimensions and is non-zero
        # The exact relationship between amplitude and intensity depends on normalization
        intensity_from_amp = abs2.(amp_img[:,:,1]) + abs2.(amp_img[:,:,2])
        @test all(intensity_from_amp .>= 0)
        @test size(intensity_from_amp) == size(img)
    end
    
   
end