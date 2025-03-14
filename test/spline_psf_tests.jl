# test/spline_psf_tests.jl

@testset "SplinePSF" begin
    @testset "2D Airy PSF Spline" begin
        # Create an Airy PSF (standard diffraction-limited PSF)
        na = 1.4
        wavelength = 0.532  # 532 nm
        airy_psf = Airy2D(na, wavelength)
        
        # Create a spline representation with 0.05μm sampling
        x_range = y_range = -1.0:0.01:1.0
        spline_psf = SplinePSF(airy_psf, x_range, y_range)
        
        # Test at various distances from center
        r_values = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8]
        
        for r in r_values
            # Test along both x and y axes
            @test isapprox(spline_psf(r, 0.0), airy_psf(r, 0.0), rtol=0.02)
            @test isapprox(spline_psf(0.0, r), airy_psf(0.0, r), rtol=0.02)
            
            # Test along diagonal
            r_diag = r / sqrt(2)
            @test isapprox(spline_psf(r_diag, r_diag), airy_psf(r_diag, r_diag), rtol=0.02)
        end
        
        # Test that spline PSF correctly represents the Airy pattern rings
        # First Airy zero is at approx 1.22 * λ / (2*NA) ≈ 0.23μm for NA=1.4, λ=0.532μm
        first_zero = 1.22 * wavelength / (2 * na)
        
        # Value at center should be high
        @test spline_psf(0.0, 0.0) > 0.9 * airy_psf(0.0, 0.0)
        
        # Value at first zero should be close to zero
        @test spline_psf(first_zero, 0.0) < 0.1
        
        # Value at second peak should be positive again
        @test spline_psf(1.5 * first_zero, 0.0) > 0.1
    end
    
    @testset "3D Scalar PSF Spline" begin
        # Create a scalar 3D PSF
        na = 1.2
        wavelength = 0.532
        n_medium = 1.33
        scalar_psf = Scalar3DPSF(na, wavelength, n_medium)
        
        # Create a spline representation with 0.05μm xy sampling and 0.1μm z sampling
        x_range = y_range = -1.0:0.05:1.0
        z_range = -1.0:0.1:1.0
        
        spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)
        
        # Test at focus (z=0)
        r_values = [0.0, 0.1, 0.3, 0.5, 0.8]
        
        for r in r_values
            # Test at focus along x-axis
            @test isapprox(spline_psf(r, 0.0, 0.0), scalar_psf(r, 0.0, 0.0), rtol=0.05)
            
            # Test at focus along diagonal
            r_diag = r / sqrt(2)
            @test isapprox(spline_psf(r_diag, r_diag, 0.0), scalar_psf(r_diag, r_diag, 0.0), rtol=0.05)
        end
        
        # Test at different z positions
        z_values = [-0.5, -0.3, 0.0, 0.3, 0.5]
        
        for z in z_values
            # Test on-axis psf
            @test isapprox(spline_psf(0.0, 0.0, z), scalar_psf(0.0, 0.0, z), rtol=0.05)
            
            # Test off-axis psf
            @test isapprox(spline_psf(0.3, 0.0, z), scalar_psf(0.3, 0.0, z), rtol=0.05)
        end
        
        # Test that the spline PSF correctly captures the PSF axial elongation
        # The PSF should be more elongated in z than in xy
        # Compare widths at which the intensity drops to 50% of maximum
        
        # Lateral profile
        xy_hm_dist = 0.0
        max_val = spline_psf(0.0, 0.0, 0.0)
        for r in 0.0:0.01:0.5
            if spline_psf(r, 0.0, 0.0) <= 0.5 * max_val
                xy_hm_dist = r
                break
            end
        end
        
        # Axial profile
        z_hm_dist = 0.0
        for z in 0.0:0.01:0.5
            if spline_psf(0.0, 0.0, z) <= 0.5 * max_val
                z_hm_dist = z
                break
            end
        end
        
        # For typical high-NA microscopes, the PSF is ~3-5x more elongated in z
        @test z_hm_dist > 2.5 * xy_hm_dist
    end
    
    @testset "Pixel Integration" begin
        # Create a 2D Airy PSF
        na = 1.4
        wavelength = 0.532
        airy_psf = Airy2D(na, wavelength)
        
        # Create spline representation with 0.05μm sampling
        x_range = y_range = -1.0:0.05:1.0
        spline_psf = SplinePSF(airy_psf, x_range, y_range)
        
        # Create camera
        pixel_edges_x = collect(0:0.1:2.0)  # 100nm pixels
        pixel_edges_y = collect(0:0.1:2.0)
        camera = IdealCamera(pixel_edges_x, pixel_edges_y)
        
        # Create emitter
        emitter = Emitter2D(1.0, 1.0, 1000.0)  # x=1μm, y=1μm, 1000 photons
        
        # Integrate pixels using both models
        spline_img = integrate_pixels(spline_psf, camera, emitter)
        airy_img = integrate_pixels(airy_psf, camera, emitter)
        
        # Compare images
        # Total photons should be similar
        @test isapprox(sum(spline_img), sum(airy_img), rtol=0.05)
        
        # Peak position should be the same
        spline_peak = findmax(spline_img)[2]
        airy_peak = findmax(airy_img)[2]
        @test spline_peak == airy_peak
        
        # Peak value should be similar
        @test isapprox(maximum(spline_img), maximum(airy_img), rtol=0.05)
    end
end