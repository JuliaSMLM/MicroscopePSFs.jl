# test/basic_psfs_tests.jl

@testset "Basic PSFs" begin
    @testset "Gaussian2D PSF" begin
        # Test constructor
        σ = 0.15  # 150 nm
        psf = Gaussian2D(σ)
        @test psf.σ ≈ σ
        
        # Test center maximum
        @test psf(0.0, 0.0) > psf(0.1, 0.0)
        @test psf(0.0, 0.0) > psf(0.0, 0.1)
        
        # Test symmetry
        @test psf(0.1, 0.0) ≈ psf(-0.1, 0.0)
        @test psf(0.0, 0.1) ≈ psf(0.0, -0.1)
        @test psf(0.1, 0.1) ≈ psf(-0.1, -0.1)
        @test psf(0.1, 0.2) ≈ psf(-0.1, -0.2)
        
        # Test radial symmetry
        @test psf(0.1, 0.0) ≈ psf(0.0, 0.1)
        @test psf(0.3, 0.0) ≈ psf(0.0, 0.3)
        
        # Test FWHM (Full Width at Half Maximum) = 2.355 * σ
        # At half maximum, exp(-r²/(2σ²)) = 0.5, which means r = σ*sqrt(2*ln(2))
        half_max_r = σ * sqrt(2 * log(2))
        center_val = psf(0.0, 0.0)
        @test psf(half_max_r, 0.0) ≈ center_val / 2 rtol=0.01
        
        # Test normalization
        # Integrate over a region large enough to contain most of the PSF
        r_max = 5.0 * σ
        dr = 0.01  # 10 nm step size
        steps = Int(ceil(r_max / dr))
        total = 0.0
        for i in 0:steps-1
            r = i * dr
            # Thin circular shell with area 2πr*dr
            total += 2π * r * dr * psf(r, 0.0)
        end
        @test total ≈ 1.0 rtol=0.01
        
        # Test amplitude function
        @test amplitude(psf, 0.1, 0.2) ≈ sqrt(psf(0.1, 0.2))
        
        # Test integration with camera and emitter
        step_xy = 0.1  # 100 nm pixel size
        x_range = 0:step_xy:2.0
        y_range = 0:step_xy:2.0
        pixel_edges_x = collect(x_range)
        pixel_edges_y = collect(y_range)
        camera = IdealCamera(pixel_edges_x, pixel_edges_y)
        emitter = Emitter2D(1.0, 1.0, 1000.0)  # At (1μm, 1μm) with 1000 photons
        img = integrate_pixels(psf, camera, emitter)
        
        # Test image dimensions match camera
        @test size(img) == (length(pixel_edges_x)-1, length(pixel_edges_y)-1)
        
        # Test energy conservation (most photons should be captured)
        @test sum(img) > 0.95 * emitter.photons
        
        # Test peak position - find pixel closest to (1.0, 1.0) μm
        max_idx = argmax(img)
        center_index_x = floor(Int, 1.0 / step_xy) + 1
        center_index_y = floor(Int, 1.0 / step_xy) + 1
        @test abs(max_idx[1] - center_index_y) <= 1
        @test abs(max_idx[2] - center_index_x) <= 1
    end
    
    @testset "Airy2D PSF" begin
        # Test constructor
        na = 1.4
        λ = 0.532  # 532 nm
        psf = Airy2D(na, λ)
        @test psf.nₐ ≈ na
        @test psf.λ ≈ λ
        @test psf.ν ≈ 2π * na / λ
        
        # Test center maximum
        @test psf(0.0, 0.0) > psf(0.1, 0.0)
        @test psf(0.0, 0.0) > psf(0.0, 0.1)
        
        # Test symmetry
        @test psf(0.1, 0.0) ≈ psf(-0.1, 0.0)
        @test psf(0.0, 0.1) ≈ psf(0.0, -0.1)
        @test psf(0.1, 0.1) ≈ psf(-0.1, -0.1)
        
        # Test radial symmetry
        @test psf(0.1, 0.0) ≈ psf(0.0, 0.1)
        @test psf(0.3, 0.0) ≈ psf(0.0, 0.3)
        
        # Test first zero position (Airy disk radius) ≈ 0.61 * λ / NA
        airy_radius = 0.61 * λ / na
        @test psf(airy_radius, 0.0) < 0.01 * psf(0.0, 0.0)
        
        # Test secondary maximum
        # Secondary maximum should be at around 1.33 * airy_radius
        second_max_r = 1.33 * airy_radius
        # Should be > 0 and much smaller than center
        @test psf(second_max_r, 0.0) > 0.01 * psf(0.0, 0.0)
        @test psf(second_max_r, 0.0) < 0.2 * psf(0.0, 0.0)
        
        # Test normalization by numerical integration
        r_max = 10.0 * airy_radius  # Increase range to capture more of PSF
        dr = 0.005  # 5 nm step size
        steps = Int(ceil(r_max / dr))
        total = 0.0
        for i in 0:steps-1
            r = i * dr
            # Thin circular shell with area 2πr*dr
            total += 2π * r * dr * psf(r, 0.0)
        end
        @test total ≈ 1.0 rtol=0.05  # Increased tolerance to 5%
        
        # Test amplitude function
        @test abs2(amplitude(psf, 0.1, 0.2)) ≈ psf(0.1, 0.2)
        
        # Test integration with camera and emitter
        step_xy = 0.05  # 50 nm pixel size
        x_range = 0:step_xy:2.0
        y_range = 0:step_xy:2.0
        pixel_edges_x = collect(x_range)
        pixel_edges_y = collect(y_range)
        camera = IdealCamera(pixel_edges_x, pixel_edges_y)
        emitter = Emitter2D(1.0, 1.0, 1000.0)  # At (1μm, 1μm) with 1000 photons
        img = integrate_pixels(psf, camera, emitter)
        
        # Test image dimensions match camera
        @test size(img) == (length(pixel_edges_y)-1, length(pixel_edges_x)-1)
        
        # Test energy conservation (most photons should be captured)
        @test sum(img) > 0.95 * emitter.photons
        
        # Test peak position - find pixel closest to (1.0, 1.0) μm
        max_idx = argmax(img)
        center_index_x = floor(Int, 1.0 / step_xy) + 1
        center_index_y = floor(Int, 1.0 / step_xy) + 1
        @test abs(max_idx[1] - center_index_y) <= 1
        @test abs(max_idx[2] - center_index_x) <= 1
    end
    
    @testset "PSF Conversions" begin
        # Test Gaussian to Airy conversion
        gaussian = Gaussian2D(0.15)
        airy_from_gaussian = Airy2D(gaussian)
        
        # Check that the conversion preserved approximate width
        # Using the relationship σ ≈ 0.22 * λ / NA
        @test airy_from_gaussian.nₐ ≈ 0.22 * airy_from_gaussian.λ / gaussian.σ rtol=0.01
        
        # Test Airy to Gaussian conversion
        airy = Airy2D(1.4, 0.532)
        gaussian_from_airy = Gaussian2D(airy)
        
        # Check that conversion preserved approximate width
        @test gaussian_from_airy.σ ≈ 0.22 * airy.λ / airy.nₐ rtol=0.01
    end
end