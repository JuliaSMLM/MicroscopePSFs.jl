# test/support_region_tests.jl

@testset "Support Region Optimization" begin
    # Setup common test parameters
    λ = 0.532  # 532 nm
    na = 1.2
    n_medium = 1.33
    
    # Camera setup with 100 nm pixels
    pixel_size = 0.1  # 100 nm
    nx, ny = 41, 41  # 41x41 pixels for a 4μm x 4μm area
    x_edges = collect(0:pixel_size:(nx*pixel_size))
    y_edges = collect(0:pixel_size:(ny*pixel_size))
    camera = IdealCamera(x_edges, y_edges)
    
    # Create emitters at different positions
    emitters = [
        Emitter2D(1.0, 1.0, 1000.0),  # Emitter at (1μm, 1μm)
        Emitter2D(3.0, 3.0, 1000.0)   # Emitter at (3μm, 3μm)
    ]
    
    # Create a 3D emitter for 3D PSF tests
    emitter_3d = Emitter3D(2.0, 2.0, 0.0, 1000.0)
    
    @testset "GaussianPSF Support Region" begin
        psf = GaussianPSF(0.15)  # σ = 150 nm
        
        # Test 1: Default - full image
        result_full = integrate_pixels(psf, camera, emitters[1])
        @test size(result_full) == (ny, nx)
        @test sum(result_full) > 0  # Should contain intensity
        
        # Test 2: Circular support region with radius 0.5μm
        result_r05 = integrate_pixels(psf, camera, emitters[1], support=0.5)
        @test size(result_r05) == (ny, nx)  # Output size should be the same
        
        # The pixels far from the emitter should be zero in the support-limited result
        far_x, far_y = floor(Int, 3.0 / pixel_size), floor(Int, 3.0 / pixel_size)
        @test result_r05[far_y, far_x] == 0.0
        
        # Pixels near the emitter should be similar to the full result
        near_x, near_y = floor(Int, 1.0 / pixel_size), floor(Int, 1.0 / pixel_size)
        @test isapprox(result_r05[near_y, near_x], result_full[near_y, near_x], rtol=1e-10)
        
        # Test 3: Explicit rectangular support region
        rect_region = (0.8, 1.2, 0.8, 1.2)  # (x_min, x_max, y_min, y_max)
        result_rect = integrate_pixels(psf, camera, emitters[1], support=rect_region)
        @test size(result_rect) == (ny, nx)
        
        # Pixels outside the rectangle should be zero
        outside_x, outside_y = floor(Int, 0.7 / pixel_size), floor(Int, 1.0 / pixel_size)
        @test result_rect[outside_y, outside_x] == 0.0
        
        # Pixels inside should match the full result
        inside_x, inside_y = floor(Int, 1.0 / pixel_size), floor(Int, 1.0 / pixel_size)
        @test isapprox(result_rect[inside_y, inside_x], result_full[inside_y, inside_x], rtol=1e-10)
    end
    
    @testset "Multiple Emitters Support Region" begin
        psf = AiryPSF(na, λ)
        
        # Test 1: Full image with multiple emitters
        result_full = integrate_pixels(psf, camera, emitters)
        @test size(result_full) == (ny, nx)
        
        # Test 2: Circular support with multiple emitters
        result_support = integrate_pixels(psf, camera, emitters, support=0.5)
        @test size(result_support) == (ny, nx)
        
        # Sum should be lower with limited support since we're cutting off PSF tails
        @test sum(result_support) < sum(result_full)
        
        # But sum should still be high (emitters contribute most energy near their centers)
        @test sum(result_support) > 0.8 * sum(result_full)
        
        # The support region is relative to each emitter, so we should see intensity
        # near both emitters (1,1) and (3,3)
        emitter1_x, emitter1_y = floor(Int, 1.0 / pixel_size), floor(Int, 1.0 / pixel_size)
        emitter2_x, emitter2_y = floor(Int, 3.0 / pixel_size), floor(Int, 3.0 / pixel_size)
        
        @test result_support[emitter1_y, emitter1_x] > 0.0
        @test result_support[emitter2_y, emitter2_x] > 0.0
        
        # Points between the emitters should be zero
        mid_x, mid_y = floor(Int, 2.0 / pixel_size), floor(Int, 2.0 / pixel_size)
        @test result_support[mid_y, mid_x] == 0.0
    end
    
    @testset "ScalarPSF Support Region" begin
        psf = ScalarPSF(na, λ, n_medium)
        
        # Test with 3D emitter
        result_full = integrate_pixels(psf, camera, emitter_3d)
        @test size(result_full) == (ny, nx)
        
        # Test with support
        result_support = integrate_pixels(psf, camera, emitter_3d, support=0.5)
        @test size(result_support) == (ny, nx)
        
        # Emitter position should have intensity
        emitter_x, emitter_y = floor(Int, 2.0 / pixel_size), floor(Int, 2.0 / pixel_size)
        @test result_support[emitter_y, emitter_x] > 0.0
        
        # Far point should be zero
        far_x, far_y = floor(Int, 0.2 / pixel_size), floor(Int, 0.2 / pixel_size)
        @test result_support[far_y, far_x] == 0.0
    end
    
    @testset "Support Region with Complex Amplitude" begin
        psf = ScalarPSF(na, λ, n_medium)
        
        # Test with amplitude instead of intensity
        amp_full = integrate_pixels_amplitude(psf, camera, emitter_3d)
        @test size(amp_full) == (ny, nx)
        @test eltype(amp_full) <: Complex
        
        # Test with support
        amp_support = integrate_pixels_amplitude(psf, camera, emitter_3d, support=0.5)
        @test size(amp_support) == (ny, nx)
        @test eltype(amp_support) <: Complex
        
        # Emitter position should have non-zero amplitude
        emitter_x, emitter_y = floor(Int, 2.0 / pixel_size), floor(Int, 2.0 / pixel_size)
        @test abs(amp_support[emitter_y, emitter_x]) > 0.0
        
        # Far point should be zero
        far_x, far_y = floor(Int, 0.2 / pixel_size), floor(Int, 0.2 / pixel_size)
        @test abs(amp_support[far_y, far_x]) == 0.0
    end
    
    @testset "SplinePSF Support Region" begin
        # Create a SplinePSF from a GaussianPSF
        gauss = GaussianPSF(0.15)
        x_range = y_range = range(-1.0, 1.0, length=41)
        spline_psf = SplinePSF(gauss, x_range, y_range)
        
        # Test with and without support
        result_full = integrate_pixels(spline_psf, camera, emitters[1])
        result_support = integrate_pixels(spline_psf, camera, emitters[1], support=0.5)
        
        @test size(result_full) == (ny, nx)
        @test size(result_support) == (ny, nx)
        
        # Emitter position should have intensity
        emitter_x, emitter_y = floor(Int, 1.0 / pixel_size), floor(Int, 1.0 / pixel_size)
        @test result_full[emitter_y, emitter_x] > 0.0
        @test result_support[emitter_y, emitter_x] > 0.0
        
        # Far point should be zero
        far_x, far_y = floor(Int, 3.0 / pixel_size), floor(Int, 3.0 / pixel_size)
        @test result_support[far_y, far_x] == 0.0
    end
    
    @testset "Subsampling Parameter" begin
        psf = GaussianPSF(0.15)
        
        # Test different sampling rates
        result_s1 = integrate_pixels(psf, camera, emitters, support=0.5, sampling=1)
        result_s2 = integrate_pixels(psf, camera, emitters, support=0.5, sampling=2)
        result_s4 = integrate_pixels(psf, camera, emitters, support=0.5, sampling=4)
        
        @test size(result_s1) == size(result_s2) == size(result_s4) == (ny, nx)
        
        # Higher sampling should give more accurate results
        # Sum should converge as sampling increases
        sum_diff_21 = abs(sum(result_s2) - sum(result_s1))
        sum_diff_42 = abs(sum(result_s4) - sum(result_s2))
        
        # The difference should decrease with higher sampling
        @test sum_diff_42 < sum_diff_21
    end
end