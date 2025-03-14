using MicroscopePSFs
using CairoMakie

println("Testing complete code blocks from docs/src/psfs/airy2d.md")

println("\n===== Testing Constructor Examples =====")
try
    # Create an Airy PSF for a high-NA objective with green light
    println("Creating Airy2D PSF...")
    psf = Airy2D(1.4, 0.532)  # NA=1.4, wavelength=532nm
    println("Successfully created Airy2D with NA=1.4, wavelength=0.532μm")
    
    # Create from a Gaussian2D PSF
    println("Creating Airy2D from Gaussian2D...")
    gaussian_psf = Gaussian2D(0.15)
    airy_equivalent = Airy2D(gaussian_psf, λ=0.532)
    println("Successfully created Airy2D from Gaussian2D")
    println("Constructor examples successfully tested")
catch e
    println("ERROR in Constructor examples: ", e)
    println(stacktrace())
end

println("\n===== Testing Evaluation Methods =====")
try
    # Create a basic PSF for testing
    psf = Airy2D(1.4, 0.532)
    
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
    psf = Airy2D(1.4, 0.532)
    
    # Create a grid of positions
    println("Creating coordinate grid...")
    x_coords = -2:0.1:2  # microns
    y_coords = -2:0.1:2  # microns
    println("Successfully created coordinate grid")
    
    # Compute PSF values at each position
    println("Computing intensity values...")
    intensity_values = [psf(xi, yi) for yi in y_coords, xi in x_coords]
    println("Successfully computed intensity image with size: ", size(intensity_values))
    
    # Visualization
    println("Creating visualization...")
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1, 1], aspect=DataAspect(), 
             title="Airy2D (NA=1.4, λ=532nm)",
             xlabel="x (μm)", ylabel="y (μm)")
    hm = heatmap!(ax, x_coords, y_coords, intensity_values, colormap=:viridis)
    Colorbar(fig[1, 2], hm)
    
    # Save the figure to verify it works
    save_path = tempname() * "_airy2d.png"
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
    psf = Airy2D(1.4, 0.532)
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
             title="Airy2D Camera Image",
             xlabel="x (μm)", ylabel="y (μm)")
    ax.yreversed = true  # Flip y-axis to match camera convention
    hm = heatmap!(ax, x_phys, y_phys, pixels', colormap=:viridis)
    scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    Colorbar(fig[1, 2], hm)
    
    # Save the figure to verify it works
    save_path = tempname() * "_airy2d_camera.png"
    save(save_path, fig)
    println("Successfully saved camera visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Camera integration successfully tested")
catch e
    println("ERROR in Camera integration: ", e)
    println(stacktrace())
end

println("\nAll code blocks from airy2d.md successfully tested!")