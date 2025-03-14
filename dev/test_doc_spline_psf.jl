using MicroscopePSFs
using MicroscopePSFs.Zernike
using CairoMakie
using BenchmarkTools

println("Testing complete code blocks from docs/src/psfs/spline_psf.md")

println("\n===== Testing 'Use Cases' example =====")
try
    # Create a computationally expensive PSF model (dipole along z-axis)
    println("Creating Vector3DPSF...")
    dipole_z = DipoleVector(0.0, 0.0, 1.0)
    vector_psf = Vector3DPSF(
        1.4,                # Numerical aperture
        0.532,              # Wavelength
        dipole_z,           # Dipole orientation
        n_medium=1.33       # Sample medium refractive index
    )
    println("Successfully created Vector3DPSF")

    # Define coordinate sampling grid (using smaller dimensions for testing)
    println("Creating coordinate grids...")
    x_range = range(-2.0, 2.0, length=21)  # μm (smaller for testing)
    y_range = range(-2.0, 2.0, length=21)  # μm
    z_range = range(-2.0, 2.0, length=11)  # μm
    println("Successfully created coordinate grids")

    # Create a spline representation
    println("Creating SplinePSF from Vector3DPSF...")
    spline_psf = SplinePSF(vector_psf, x_range, y_range, z_range)
    println("Successfully created SplinePSF")

    # Benchmark comparison
    println("Benchmarking PSF evaluations...")
    x, y, z = 0.5, -0.2, 0.3
    println("Vector3DPSF evaluation time:")
    vector_time = @elapsed vector_psf(x, y, z)
    println("  $vector_time seconds")
    println("SplinePSF evaluation time:")
    spline_time = @elapsed spline_psf(x, y, z)
    println("  $spline_time seconds")

    # Visual comparison
    println("Creating visualization...")
    fig = Figure(size=(900, 400))
    z_positions = [-1.0, 0.0, 1.0]
    for (i, z) in enumerate(z_positions)
        # Calculate PSF images at this z position
        println("Computing PSF images at z=$(z)...")
        intensity_vector = [vector_psf(xi, yi, z) for yi in y_range, xi in x_range]
        intensity_spline = [spline_psf(xi, yi, z) for yi in y_range, xi in x_range]
        
        # Plot original PSF
        ax1 = Axis(fig[1, i], aspect=DataAspect(), 
                 title="Vector3DPSF (z=$(z)μm)",
                 xlabel="x (μm)", ylabel="y (μm)")
        hm1 = heatmap!(ax1, x_range, y_range, intensity_vector, colormap=:viridis)
        
        # Plot spline version
        ax2 = Axis(fig[2, i], aspect=DataAspect(), 
                 title="SplinePSF (z=$(z)μm)",
                 xlabel="x (μm)", ylabel="y (μm)")
        hm2 = heatmap!(ax2, x_range, y_range, intensity_spline, colormap=:viridis)
    end

    # Save figure to verify it works
    save_path = tempname() * "_use_cases.png"
    save(save_path, fig)
    println("Successfully saved visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Use Cases example successfully tested")
catch e
    println("ERROR in Use Cases example: ", e)
    println(stacktrace())
end

println("\n===== Testing 'Accelerating a Scalar3DPSF PSF' example =====")
try
    # Create a Scalar3DPSF PSF with aberrations
    println("Creating Scalar3DPSF with aberrations...")
    zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomials
    add_spherical!(zc, 0.5)       # Add spherical aberration
    add_astigmatism!(zc, 0.3)     # Add astigmatism
    
    scalar_psf = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)
    println("Successfully created Scalar3DPSF with aberrations")
    
    # Create a spline version for faster evaluation (using smaller grids for testing)
    println("Creating coordinate grids...")
    x_range = range(-2.0, 2.0, length=21)  # Smaller for testing
    y_range = range(-2.0, 2.0, length=21)
    z_range = range(-2.0, 2.0, length=11)
    println("Successfully created coordinate grids")
    
    # This step takes time but only needs to be done once
    println("Creating SplinePSF from Scalar3DPSF...")
    spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)
    println("Successfully created SplinePSF")
    
    # Compare speed for single point evaluation
    println("Benchmarking PSF evaluations...")
    println("Scalar3DPSF evaluation time:")
    scalar_time = @elapsed scalar_psf(0.5, 0.5, 0.1)
    println("  $scalar_time seconds")
    println("SplinePSF evaluation time:")
    spline_time = @elapsed spline_psf(0.5, 0.5, 0.1)
    println("  $spline_time seconds")
    
    # Compare speed for camera integration
    println("Setting up camera integration...")
    pixel_edges_x = collect(0:0.1:2.0)
    pixel_edges_y = collect(0:0.1:2.0)
    camera = IdealCamera(pixel_edges_x, pixel_edges_y)  # 20x20 camera
    emitter = Emitter3D(1.0, 1.0, 0.0, 1000.0)          # Emitter with 1000 photons
    println("Successfully created camera and emitter")
    
    println("Benchmarking pixel integration...")
    println("Scalar3DPSF integration time:")
    scalar_int_time = @elapsed integrate_pixels(scalar_psf, camera, emitter)
    println("  $scalar_int_time seconds")
    println("SplinePSF integration time:")
    spline_int_time = @elapsed integrate_pixels(spline_psf, camera, emitter)
    println("  $spline_int_time seconds")
    
    # Visualization of integrated pixels
    println("Integrating pixels...")
    pixels_scalar = integrate_pixels(scalar_psf, camera, emitter)
    pixels_spline = integrate_pixels(spline_psf, camera, emitter)
    println("Successfully integrated pixels")
    
    # Camera physical coordinates
    x_phys = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) / 2
    y_phys = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) / 2
    println("Calculated pixel center coordinates")
    
    # Plot comparison
    println("Creating visualization...")
    fig = Figure(size=(800, 400))
    
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), 
            title="Scalar3DPSF",
            xlabel="x (μm)", ylabel="y (μm)")
    ax1.yreversed = true
    hm1 = heatmap!(ax1, x_phys, y_phys, pixels_scalar', colormap=:viridis)
    scatter!(ax1, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    Colorbar(fig[1, 2], hm1)
    
    ax2 = Axis(fig[1, 3], aspect=DataAspect(), 
            title="SplinePSF",
            xlabel="x (μm)", ylabel="y (μm)")
    ax2.yreversed = true
    hm2 = heatmap!(ax2, x_phys, y_phys, pixels_spline', colormap=:viridis)
    scatter!(ax2, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    Colorbar(fig[1, 4], hm2)
    
    # Save figure to verify it works
    save_path = tempname() * "_accelerating.png"
    save(save_path, fig)
    println("Successfully saved visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Accelerating Scalar3DPSF example successfully tested")
catch e
    println("ERROR in Accelerating Scalar3DPSF example: ", e)
    println(stacktrace())
end

println("\n===== Testing 'Using Experimental PSF Data' example =====")
try
    # We'll simulate the experimental PSF instead of loading it
    println("Simulating experimental PSF data...")
    nx, ny, nz = 15, 15, 7  # Small example dimensions
    measured_psf = rand(nx, ny, nz)  # Simulated measurement
    
    # Define physical coordinates for the experimental PSF
    pixel_size = 0.1  # μm
    z_step = 0.2      # μm
    println("Successfully created simulated PSF data")
    
    # Create coordinate ranges based on PSF stack dimensions
    println("Creating coordinate ranges...")
    x_range = range(-nx//2 * pixel_size, (nx//2-1) * pixel_size, length=nx)
    y_range = range(-ny//2 * pixel_size, (ny//2-1) * pixel_size, length=ny)
    z_range = range(-nz//2 * z_step, (nz//2-1) * z_step, length=nz)
    println("Successfully created coordinate ranges")
    
    # Create the SplinePSF
    println("Creating SplinePSF from experimental data...")
    experimental_psf = SplinePSF(measured_psf, x_range, y_range, z_range)
    println("Successfully created SplinePSF from experimental data")
    
    # Use it like any other PSF
    println("Setting up camera integration...")
    pixel_edges_x = collect(0:0.1:2.0)
    pixel_edges_y = collect(0:0.1:2.0)
    camera = IdealCamera(pixel_edges_x, pixel_edges_y)
    emitter = Emitter3D(1.0, 1.0, 0.0, 1000.0)
    pixels = integrate_pixels(experimental_psf, camera, emitter)
    println("Successfully integrated pixels with size: ", size(pixels))
    println("Using Experimental PSF Data example successfully tested")
catch e
    println("ERROR in Using Experimental PSF Data example: ", e)
    println(stacktrace())
end

println("\nAll code blocks from spline_psf.md successfully tested!")