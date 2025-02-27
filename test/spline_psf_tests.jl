# test/spline_psf_tests.jl

@testset "SplinePSF" begin
    # Test creation from data
    @testset "Creation from data" begin
        # Create a simple Gaussian pattern
        sz = 21
        x = y = z = range(-1.0, 1.0, length=sz)
        sigma = 0.3
        
        # Generate 3D Gaussian
        data = zeros(sz, sz, sz)
        for iz in 1:sz, ix in 1:sz, iy in 1:sz
            r2 = x[ix]^2 + y[iy]^2 + z[iz]^2
            data[iy, ix, iz] = exp(-r2 / (2*sigma^2))
        end
        
        # Create SplinePSF
        spline_psf = SplinePSF(data; x_coords=x, y_coords=y, z_coords=z)
        
        # Test interpolation accuracy
        interp_point = (0.123, -0.456, 0.789)
        exact_value = exp(-(interp_point[1]^2 + interp_point[2]^2 + interp_point[3]^2) / (2*sigma^2))
        interp_value = spline_psf(interp_point...)
        
        @test isapprox(interp_value * sum(data), exact_value, rtol=1e-3)
    end
    
    # Test creation from analytical PSF
    @testset "Creation from analytical PSF" begin
        # Create Gaussian PSF
        gauss_psf = Gaussian2D(0.15)
        
        # Convert to SplinePSF with dense sampling
        spline_psf = SplinePSF(gauss_psf; 
                              x_range=range(-0.5, 0.5, length=41),
                              y_range=range(-0.5, 0.5, length=41),
                              z_range=range(-0.5, 0.5, length=5))
        
        # Test at some random positions
        for _ in 1:10
            x = rand() * 0.4 - 0.2
            y = rand() * 0.4 - 0.2
            
            # Values should be very close
            @test isapprox(spline_psf(x, y, 0.0), gauss_psf(x, y), rtol=1e-3)
        end
    end
    
    # Test 2D evaluation
    @testset "2D evaluation" begin
        # Create simple SplinePSF
        sz = 11
        x = y = range(-1.0, 1.0, length=sz)
        z = range(-0.5, 0.5, length=5)
        
        data = zeros(sz, sz, 5)
        for iz in 1:5, ix in 1:sz, iy in 1:sz
            r2 = x[ix]^2 + y[iy]^2
            data[iy, ix, iz] = exp(-r2 / 0.2)
        end
        
        spline_psf = SplinePSF(data; x_coords=x, y_coords=y, z_coords=z)
        
        # 2D evaluation should use z=0
        @test isapprox(spline_psf(0.3, 0.3), spline_psf(0.3, 0.3, 0.0), rtol=1e-10)
    end
    
    # Test amplitude
    @testset "Amplitude" begin
        # Create simple SplinePSF
        sz = 11
        x = y = range(-1.0, 1.0, length=sz)
        z = range(-0.5, 0.5, length=3)
        
        data = zeros(sz, sz, 3)
        for iz in 1:3, ix in 1:sz, iy in 1:sz
            r2 = x[ix]^2 + y[iy]^2
            data[iy, ix, iz] = exp(-r2 / 0.2)
        end
        
        spline_psf = SplinePSF(data; x_coords=x, y_coords=y, z_coords=z)
        
        # Amplitude should be sqrt of intensity
        @test isapprox(abs2(amplitude(spline_psf, 0.3, 0.4, 0.1)), spline_psf(0.3, 0.4, 0.1), rtol=1e-10)
    end
    
    # Test pixel integration
    @testset "Pixel integration" begin
        # Create simple SplinePSF
        sz = 31
        x = y = range(-1.0, 1.0, length=sz)
        z = [0.0]
        
        data = zeros(sz, sz, 1)
        for ix in 1:sz, iy in 1:sz
            r2 = x[ix]^2 + y[iy]^2
            data[iy, ix, 1] = exp(-r2 / 0.2)
        end
        
        spline_psf = SplinePSF(data; x_coords=x, y_coords=y, z_coords=z)
        
        # Create camera
        camera = IdealCamera(0:0.1:2.0, 0:0.1:2.0)  # 20x20 pixels, 100nm each
        
        # Create emitter
        emitter = Emitter2D(1.0, 1.0, 1000.0)  # x=1μm, y=1μm, 1000 photons
        
        # Integrate pixels
        pixels = integrate_pixels(spline_psf, camera, emitter)
        
        # Check that result normalizes to 1
        @test isapprox(sum(pixels), 1.0, rtol=1e-10)
        
        # Highest intensity should be near emitter position
        max_idx = findmax(pixels)[2]
        ctr_x, ctr_y = 10, 10  # Pixel indices near (1,1)
        @test abs(max_idx[1] - ctr_y) <= 1 && abs(max_idx[2] - ctr_x) <= 1
    end
    
    # Test save/load
    @testset "Save and load" begin
        # Create a simple SplinePSF
        sz = 11
        x = y = range(-1.0, 1.0, length=sz)
        z = range(-0.5, 0.5, length=3)
        
        data = zeros(sz, sz, 3)
        for iz in 1:3, ix in 1:sz, iy in 1:sz
            r2 = x[ix]^2 + y[iy]^2 + z[iz]^2
            data[iy, ix, iz] = exp(-r2 / 0.2)
        end
        
        original_psf = SplinePSF(data; x_coords=x, y_coords=y, z_coords=z)
        
        # Save to temporary file
        temp_file = tempname() * ".h5"
        save_spline_psf(temp_file, original_psf)
        
        # Load from file
        loaded_psf = load_spline_psf(temp_file)
        
        # Compare values at random points
        for _ in 1:10
            test_point = (rand()*1.8-0.9, rand()*1.8-0.9, rand()*0.9-0.45)
            @test isapprox(original_psf(test_point...), loaded_psf(test_point...), rtol=1e-3)
        end
        
        # Clean up
        rm(temp_file)
    end
    
    # Test automatic differentiation with Zygote
    @testset "Automatic differentiation" begin
        # Create a simple SplinePSF
        sz = 21
        x = y = range(-1.0, 1.0, length=sz)
        z = [0.0] # 2D for simplicity
        
        data = zeros(sz, sz, 1)
        sigma = 0.3
        for ix in 1:sz, iy in 1:sz
            r2 = x[ix]^2 + y[iy]^2
            data[iy, ix, 1] = exp(-r2 / (2*sigma^2))
        end
        
        spline_psf = SplinePSF(data; x_coords=x, y_coords=y, z_coords=z)
        
        # Define a function to differentiate
        f(p) = spline_psf(p[1], p[2], 0.0)
        
        # Use Zygote to get gradient
        p0 = [0.2, 0.3]
        grad = Zygote.gradient(f, p0)[1]
        
        # Compare to analytical gradient of Gaussian
        r2 = p0[1]^2 + p0[2]^2
        exact_val = exp(-r2 / (2*sigma^2))
        exact_grad = -exact_val .* p0 ./ sigma^2
        
        # Normalize for comparison (since our SplinePSF is normalized)
        norm_factor = sum(data)
        exact_grad = exact_grad / norm_factor
        
        # Test gradient accuracy
        @test isapprox(grad[1], exact_grad[1], rtol=1e-2)
        @test isapprox(grad[2], exact_grad[2], rtol=1e-2)
    end
end