using MicroscopePSFs
using CairoMakie

println("Testing code blocks from docs/src/index.md")

# Quick Start Example
println("\nTesting Quick Start example...")

# Create a Gaussian2D PSF
psf = Gaussian2D(0.15)  # sigma = 150nm
println("Created Gaussian2D PSF with σ = 0.15")

# Evaluate at a specific position
intensity = psf(0.1, 0.2)  # Intensity at (x,y) = (0.1μm, 0.2μm)
println("Intensity at (0.1, 0.2): ", intensity)

# Generate a PSF image
x = range(-1, 1, length=101)  # μm
y = range(-1, 1, length=101)  # μm
img = [psf(xi, yi) for yi in y, xi in x]
println("Generated PSF image with size: ", size(img))

# Testing Visualization
println("\nTesting PSF visualization...")
try
    # Visualize the PSF
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1, 1], aspect=DataAspect(),
            title="Gaussian PSF (σ=150nm)",
            xlabel="x (μm)", 
            ylabel="y (μm)")
    hm = heatmap!(ax, x, y, img, colormap=:viridis)
    Colorbar(fig[1, 2], hm)
    println("Successfully created PSF visualization")
    
    # Testing camera image generation
    println("\nTesting camera integration and visualization...")
    # The IdealCamera constructor should use vectors, not ranges directly
    pixel_edges_x = collect(0:0.1:2.0)  # Convert to Vector
    pixel_edges_y = collect(0:0.1:2.0)  # Convert to Vector
    camera = IdealCamera(pixel_edges_x, pixel_edges_y)  # 20x20 pixels, 100nm size
    emitter = Emitter2D(1.0, 1.0, 1000.0)               # At (1μm, 1μm) with 1000 photons
    pixels = integrate_pixels(psf, camera, emitter)
    println("Successfully created camera and integrated pixels with size: ", size(pixels))
    
    # Camera physical coordinates
    x_phys = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) / 2
    y_phys = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) / 2
    println("Calculated pixel center coordinates")
    
    # Visualize the camera image
    ax2 = Axis(fig[2, 1:2], aspect=DataAspect(),
            title="Integrated Camera Image",
            xlabel="x (μm)", 
            ylabel="y (μm)")
    ax2.yreversed = true  # Flip y axis to match camera convention
    hm2 = heatmap!(ax2, x_phys, y_phys, pixels', colormap=:viridis)
    scatter!(ax2, [emitter.x], [emitter.y], 
            color=:red, marker=:cross, markersize=15)
    println("Successfully created camera image visualization")
    
    # Save the figure to verify it worked 
    # (saving to a temporary location that will be cleaned up after testing)
    save_path = tempname() * ".png"
    save(save_path, fig)
    println("Successfully saved visualization to: ", save_path)
    rm(save_path)  # Clean up the temporary file

    # Additional test to verify we can display the figure
    # This won't actually show the figure during testing, but will check the code runs
    println("Testing display(fig)...")
    display(fig)
    println("Display command executed successfully")
    
catch e
    println("ERROR in visualization: ", e)
    println(stacktrace())
end

# Testing with direct ranges to confirm they don't work with IdealCamera
println("\nTesting IdealCamera with direct ranges (should fail)...")
try
    camera_bad = IdealCamera(0:0.1:2.0, 0:0.1:2.0)
    println("WARNING: IdealCamera constructor accepted ranges directly - this is inconsistent with our fixes")
catch e
    println("Expected error with direct ranges: ", e)
end

println("\nAll code blocks tested from index.md")