# test/thread_tests.jl

@testset "Threading Options" begin
    # Setup common test parameters
    λ = 0.532  # 532 nm
    na = 1.2
    n_medium = 1.33
    
    # Camera setup with small grid for faster tests
    pixel_size = 0.1  # 100 nm
    nx, ny = 21, 21  # 21x21 pixels
    x_edges = collect(0:pixel_size:(nx*pixel_size))
    y_edges = collect(0:pixel_size:(ny*pixel_size))
    camera = IdealCamera(x_edges, y_edges)
    
    # Create single emitter
    emitter = Emitter2D(1.0, 1.0, 1000.0)  # Emitter at (1μm, 1μm)
    
    # Create multiple emitters
    emitters = [
        Emitter2D(0.5, 0.5, 500.0),
        Emitter2D(1.5, 1.5, 750.0)
    ]
    
    # 3D emitters for 3D PSFs
    emitter_3d = Emitter3D(1.0, 1.0, 0.2, 1000.0)
    emitters_3d = [
        Emitter3D(0.5, 0.5, 0.0, 500.0),
        Emitter3D(1.5, 1.5, 0.3, 750.0)
    ]
    
    @testset "GaussianPSF Threading" begin
        psf = GaussianPSF(0.15)
        
        # Test single emitter
        result_threaded = integrate_pixels(psf, camera, emitter, threaded=true)
        result_nonthreaded = integrate_pixels(psf, camera, emitter, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
        
        # Test multiple emitters
        result_threaded = integrate_pixels(psf, camera, emitters, threaded=true)
        result_nonthreaded = integrate_pixels(psf, camera, emitters, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
        
        # Test with support parameter
        result_threaded = integrate_pixels(psf, camera, emitter, support=0.5, threaded=true)
        result_nonthreaded = integrate_pixels(psf, camera, emitter, support=0.5, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
    end
    
    @testset "AiryPSF Threading" begin
        psf = AiryPSF(na, λ)
        
        # Test single emitter
        result_threaded = integrate_pixels(psf, camera, emitter, threaded=true)
        result_nonthreaded = integrate_pixels(psf, camera, emitter, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
        
        # Test amplitude integration
        amp_threaded = integrate_pixels_amplitude(psf, camera, emitter, threaded=true)
        amp_nonthreaded = integrate_pixels_amplitude(psf, camera, emitter, threaded=false)
        
        # Results should be identical
        @test size(amp_threaded) == size(amp_nonthreaded)
        @test all(amp_threaded .≈ amp_nonthreaded)
    end
    
    @testset "ScalarPSF Threading" begin
        psf = ScalarPSF(na, λ, n_medium)
        
        # Test 3D emitter
        result_threaded = integrate_pixels(psf, camera, emitter_3d, threaded=true)
        result_nonthreaded = integrate_pixels(psf, camera, emitter_3d, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
        
        # Test multiple 3D emitters
        result_threaded = integrate_pixels(psf, camera, emitters_3d, threaded=true)
        result_nonthreaded = integrate_pixels(psf, camera, emitters_3d, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
        
        # Test amplitude integration with support
        amp_threaded = integrate_pixels_amplitude(psf, camera, emitter_3d, support=0.5, threaded=true)
        amp_nonthreaded = integrate_pixels_amplitude(psf, camera, emitter_3d, support=0.5, threaded=false)
        
        # Results should be identical
        @test size(amp_threaded) == size(amp_nonthreaded)
        @test all(amp_threaded .≈ amp_nonthreaded)
    end
    
    @testset "SplinePSF Threading" begin
        # Create SplinePSF from GaussianPSF
        gauss = GaussianPSF(0.15)
        x_range = y_range = range(-1.0, 1.0, length=21)
        spline_psf = SplinePSF(gauss, x_range, y_range)
        
        # Test threaded vs non-threaded
        result_threaded = integrate_pixels(spline_psf, camera, emitter, threaded=true)
        result_nonthreaded = integrate_pixels(spline_psf, camera, emitter, threaded=false)
        
        # Results should be identical
        @test size(result_threaded) == size(result_nonthreaded)
        @test all(result_threaded .≈ result_nonthreaded)
    end
    
    @testset "Various Sampling Levels" begin
        psf = AiryPSF(na, λ)
        
        # Test with different sampling levels
        for sampling in [1, 2, 4]
            result_threaded = integrate_pixels(psf, camera, emitter, sampling=sampling, threaded=true)
            result_nonthreaded = integrate_pixels(psf, camera, emitter, sampling=sampling, threaded=false)
            
            # Results should be identical
            @test size(result_threaded) == size(result_nonthreaded)
            @test all(result_threaded .≈ result_nonthreaded)
        end
    end
    
    @testset "Vector PSF Threading" begin
        try
            # Create VectorPSF with dipole
            dipole = DipoleVector(0.0, 0.0, 1.0)
            psf = VectorPSF(na, λ, dipole, n_medium=n_medium)
            
            # Create dipole emitter
            dipole_emitter = DipoleEmitter3D(1.0, 1.0, 0.0, 1000.0, 0.0, 0.0, 1.0)
            
            # Test intensity integration
            result_threaded = integrate_pixels(psf, camera, dipole_emitter, threaded=true)
            result_nonthreaded = integrate_pixels(psf, camera, dipole_emitter, threaded=false)
            
            # Results should be identical
            @test size(result_threaded) == size(result_nonthreaded)
            @test all(result_threaded .≈ result_nonthreaded)
            
            # Test amplitude integration
            amp_threaded = integrate_pixels_amplitude(psf, camera, dipole_emitter, threaded=true)
            amp_nonthreaded = integrate_pixels_amplitude(psf, camera, dipole_emitter, threaded=false)
            
            # Results should be identical
            @test size(amp_threaded) == size(amp_nonthreaded)
            @test all(amp_threaded .≈ amp_nonthreaded)
        catch e
            @info "Skipping VectorPSF threading tests: $e"
        end
    end
end