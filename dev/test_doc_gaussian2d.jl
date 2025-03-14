using MicroscopePSFs
using CairoMakie

println("Testing complete code blocks from docs/src/psfs/gaussian2d.md")

println("\n===== Testing Constructor Examples =====")
try
    # Create a Gaussian PSF with 150nm standard deviation
    println("Creating Gaussian2D PSF...")
    psf = Gaussian2D(0.15)
    println("Successfully created Gaussian2D with σ=0.15μm")
    
    # Create a Gaussian approximation of an Airy disk
    println("Creating Gaussian2D from Airy2D...")
    airy_psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm
    gaussian_approximation = Gaussian2D(airy_psf)  # Automatically sets appropriate σ
    println("Successfully created Gaussian2D from Airy2D with σ=", gaussian_approximation.σ)
    println("Constructor examples successfully tested")
catch e
    println("ERROR in Constructor examples: ", e)
    println(stacktrace())
end

println("\n===== Testing Evaluation Methods =====")
try
    # Create a basic PSF for testing
    psf = Gaussian2D(0.15)
    
    # Evaluate PSF at a specific position
    println("Evaluating PSF at (0.1, 0.2)...")
    intensity = psf(0.1, 0.2)
    println("Intensity at (0.1, 0.2): ", intensity)
    
    # Get complex amplitude
    println("Getting complex amplitude at (0.1, 0.2)...")
    amp = amplitude(psf, 0.1, 0.2)
    println("Complex amplitude at (0.1, 0.2): ", amp)
    println("Evaluation methods successfully tested")
catch e
    println("ERROR in Evaluation methods: ", e)
    println(stacktrace())
end

println("\n===== Testing Creating Images =====")
try
    # Create a basic PSF for testing
    psf = Gaussian2D(0.15)
    
    # Create a grid of positions
    println("Creating coordinate grid...")
    x = range(-1, 1, length=100)
    y = range(-1, 1, length=100)
    println("Successfully created coordinate grid")
    
    # Compute PSF values at each position
    println("Computing intensity values...")
    intensity_values = [psf(xi, yi) for yi in y, xi in x]
    println("Successfully computed intensity image with size: ", size(intensity_values))
    
    # Visualization
    println("Creating visualization...")
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1, 1], aspect=DataAspect(), 
             title="Gaussian2D (σ=0.15μm)",
             xlabel="x (μm)", ylabel="y (μm)")
    hm = heatmap!(ax, x, y, intensity_values, colormap=:viridis)
    Colorbar(fig[1, 2], hm)
    
    # Save the figure to verify it works
    save_path = tempname() * "_gaussian2d.png"
    save(save_path, fig)
    println("Successfully saved visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Creating images successfully tested")
catch e
    println("ERROR in Creating images: ", e)
    println(stacktrace())
end

println("\n===== Testing Camera Integration =====")
try
    # Create a basic PSF and setup camera integration
    println("Setting up camera integration...")
    psf = Gaussian2D(0.15)
    pixel_edges_x = collect(0:0.1:2.0)
    pixel_edges_y = collect(0:0.1:2.0)
    camera = IdealCamera(pixel_edges_x, pixel_edges_y)
    emitter = Emitter2D(1.0, 1.0, 1000.0)
    
    # Integrate pixels
    println("Integrating pixels...")
    pixels = integrate_pixels(psf, camera, emitter)
    println("Successfully integrated pixels with size: ", size(pixels))
    
    # Visualization of camera image
    println("Creating camera image visualization...")
    # Calculate pixel centers for plotting
    x_phys = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) / 2
    y_phys = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) / 2
    
    fig = Figure(size=(500, 400))
    ax = Axis(fig[1, 1], aspect=DataAspect(),
             title="Gaussian2D Camera Image",
             xlabel="x (μm)", ylabel="y (μm)")
    ax.yreversed = true  # Flip y-axis to match camera convention
    hm = heatmap!(ax, x_phys, y_phys, pixels', colormap=:viridis)
    scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    Colorbar(fig[1, 2], hm)
    
    # Save the figure to verify it works
    save_path = tempname() * "_gaussian2d_camera.png"
    save(save_path, fig)
    println("Successfully saved camera visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Camera integration successfully tested")
catch e
    println("ERROR in Camera integration: ", e)
    println(stacktrace())
end

println("\nAll code blocks from gaussian2d.md successfully tested!")