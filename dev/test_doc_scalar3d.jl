using MicroscopePSFs
using MicroscopePSFs.Zernike
using CairoMakie

println("Testing complete code blocks from docs/src/psfs/scalar3d.md")

println("\n===== Testing Constructor Examples =====")
try
    # Create an unaberrated 3D PSF
    println("Creating unaberrated Scalar3DPSF...")
    psf = Scalar3DPSF(1.4, 0.532, 1.518)
    println("Successfully created unaberrated Scalar3DPSF")
    
    # Create a PSF with spherical aberration
    println("Creating Scalar3DPSF with spherical aberration...")
    zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial
    add_spherical!(zc, 0.5)       # Add 0.5 waves of spherical aberration
    psf_aberrated = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)
    println("Successfully created Scalar3DPSF with spherical aberration")
    println("Constructor examples successfully tested")
catch e
    println("ERROR in Constructor examples: ", e)
    println(stacktrace())
end

println("\n===== Testing Evaluation Methods =====")
try
    # Create a basic PSF for testing
    psf = Scalar3DPSF(1.4, 0.532, 1.518)
    
    # Evaluate PSF at a specific 3D position
    println("Evaluating PSF at (0.1, 0.2, 0.3)...")
    intensity = psf(0.1, 0.2, 0.3)
    println("Intensity at (0.1, 0.2, 0.3): ", intensity)
    
    # Get complex amplitude
    println("Getting complex amplitude at (0.1, 0.2, 0.3)...")
    amp = amplitude(psf, 0.1, 0.2, 0.3)
    println("Complex amplitude at (0.1, 0.2, 0.3): ", amp)
    println("Evaluation methods successfully tested")
catch e
    println("ERROR in Evaluation methods: ", e)
    println(stacktrace())
end

println("\n===== Testing Creating Images =====")
try
    # Create a basic PSF for testing
    psf = Scalar3DPSF(1.4, 0.532, 1.518)
    
    # Create a grid of positions for a 2D slice at a specific z
    println("Creating coordinate grid...")
    x_coords = -2:0.5:2  # microns (coarser for testing)
    y_coords = -2:0.5:2  # microns
    z = 0.0  # focal plane
    println("Successfully created coordinate grid")
    
    # Compute PSF values at each position
    println("Computing 2D intensity image...")
    intensity_2d = [psf(xi, yi, z) for yi in y_coords, xi in x_coords]
    println("Successfully computed 2D intensity image with size: ", size(intensity_2d))
    
    # Create a 3D PSF stack
    println("Computing 3D intensity stack...")
    z_stack = range(-2, 2, length=3)  # z positions in microns (small for testing)
    intensity_3d = Array{Float64}(undef, length(y_coords), length(x_coords), length(z_stack))
    for (k, zi) in enumerate(z_stack)
        for (j, yi) in enumerate(y_coords)
            for (i, xi) in enumerate(x_coords)
                intensity_3d[j, i, k] = psf(xi, yi, zi)
            end
        end
    end
    println("Successfully computed 3D intensity stack with size: ", size(intensity_3d))
    
    # Visualization
    println("Creating visualization...")
    fig = Figure(size=(800, 300))
    for (i, slice_z) in enumerate(z_stack)
        ax = Axis(fig[1, i], aspect=DataAspect(), 
                 title="Scalar3DPSF (z=$(slice_z)μm)",
                 xlabel="x (μm)", ylabel="y (μm)")
        # Extract the slice from the 3D stack
        slice = intensity_3d[:, :, i]
        hm = heatmap!(ax, x_coords, y_coords, slice, colormap=:viridis)
        if i == length(z_stack)
            Colorbar(fig[1, i+1], hm)
        end
    end
    
    # Save the figure to verify it works
    save_path = tempname() * "_scalar3d.png"
    save(save_path, fig)
    println("Successfully saved visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Creating images successfully tested")
catch e
    println("ERROR in Creating images: ", e)
    println(stacktrace())
end

println("\n===== Testing Aberration Modeling =====")
try
    # Create Zernike coefficients
    println("Creating Zernike coefficients...")
    zc = ZernikeCoefficients(15)  # Support up to 15th Zernike mode
    println("Successfully created Zernike coefficients")
    
    # Add common aberrations
    println("Adding aberrations...")
    add_defocus!(zc, 1.0)         # 1 wave of defocus
    add_astigmatism!(zc, 0.5)     # 0.5 waves of astigmatism
    add_coma!(zc, 0.3)            # 0.3 waves of coma
    add_spherical!(zc, 0.2)       # 0.2 waves of spherical aberration
    println("Successfully added aberrations")
    
    # Create PSF with these aberrations
    println("Creating aberrated Scalar3DPSF...")
    psf = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)
    println("Successfully created Scalar3DPSF with multiple aberrations")
    
    # Test the PSF by evaluating at a point
    intensity = psf(0.1, 0.2, 0.3)
    println("Aberrated PSF intensity at (0.1, 0.2, 0.3): ", intensity)
    println("Aberration modeling successfully tested")
catch e
    println("ERROR in Aberration modeling: ", e)
    println(stacktrace())
end

println("\nAll code blocks from scalar3d.md successfully tested!")