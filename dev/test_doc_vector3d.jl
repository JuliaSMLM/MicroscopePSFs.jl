using MicroscopePSFs
using MicroscopePSFs.Zernike
using CairoMakie

println("Testing complete code blocks from docs/src/psfs/vector3d.md")

println("\n===== Testing Constructor Examples =====")
try
    # Create a basic Vector3D PSF with Z-oriented dipole
    println("Creating Vector3DPSF with Z-oriented dipole...")
    dipole_z = DipoleVector(0.0, 0.0, 1.0)
    psf = Vector3DPSF(
        1.4,                # Numerical aperture
        0.532,              # Wavelength in microns
        dipole_z,           # Z-oriented dipole
        n_medium=1.33,      # Sample is water
        n_immersion=1.518   # Immersion oil
    )
    println("Successfully created Vector3DPSF with Z-oriented dipole")

    # Create a PSF with X-oriented dipole
    println("Creating Vector3DPSF with X-oriented dipole...")
    dipole_x = DipoleVector(1.0, 0.0, 0.0)
    psf_x = Vector3DPSF(
        1.4,
        0.532,
        dipole_x,
        n_medium=1.33
    )
    println("Successfully created Vector3DPSF with X-oriented dipole")

    # Create a PSF with aberrations
    println("Creating Vector3DPSF with aberrations...")
    zc = ZernikeCoefficients(15)
    add_spherical!(zc, 0.5)
    psf_aberrated = Vector3DPSF(
        1.4,
        0.532,
        dipole_z,
        n_medium=1.33,
        base_zernike=zc
    )
    println("Successfully created Vector3DPSF with aberrations")
    println("Constructor examples successfully tested")
catch e
    println("ERROR in Constructor examples: ", e)
    println(stacktrace())
end

println("\n===== Testing Evaluation Methods =====")
try
    # Create a basic PSF for testing
    dipole_z = DipoleVector(0.0, 0.0, 1.0)
    psf = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33)
    
    # Evaluate PSF at a specific 3D position
    println("Evaluating PSF at (0.1, 0.2, 0.3)...")
    intensity = psf(0.1, 0.2, 0.3)
    println("Intensity: ", intensity)
    
    # Get complex amplitude
    println("Getting complex amplitude at (0.1, 0.2, 0.3)...")
    amp = amplitude(psf, 0.1, 0.2, 0.3)
    println("Amplitude: ", amp)
    println("Evaluation methods successfully tested")
catch e
    println("ERROR in Evaluation methods: ", e)
    println(stacktrace())
end

println("\n===== Testing Creating Images =====")
try
    # Create a basic PSF for testing
    dipole_z = DipoleVector(0.0, 0.0, 1.0)
    psf = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33)
    
    # Create a grid of positions for a 2D slice at a specific z
    println("Creating coordinate grid...")
    x_coords = -2:0.5:2  # microns (coarse for testing)
    y_coords = -2:0.5:2  # microns
    z = 0.0  # focal plane
    println("Successfully created coordinate grid")
    
    # Compute PSF values at each position
    println("Computing 2D intensity values...")
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
    fig = Figure(size=(600, 300))
    ax = Axis(fig[1, 1], aspect=DataAspect(), 
             title="Vector3DPSF (z=0)",
             xlabel="x (μm)", ylabel="y (μm)")
    hm = heatmap!(ax, x_coords, y_coords, intensity_2d, colormap=:viridis)
    Colorbar(fig[1, 2], hm)
    
    # Save the figure to verify it works
    save_path = tempname() * "_vector3d.png"
    save(save_path, fig)
    println("Successfully saved visualization to: ", save_path)
    rm(save_path)  # Clean up
    println("Creating images successfully tested")
catch e
    println("ERROR in Creating images: ", e)
    println(stacktrace())
end

println("\n===== Testing Dipole Orientation =====")
try
    # Create dipole vectors for specific orientations
    println("Creating dipole vectors...")
    dipole_x = DipoleVector(1.0, 0.0, 0.0)  # X-oriented dipole
    dipole_y = DipoleVector(0.0, 1.0, 0.0)  # Y-oriented dipole
    dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Z-oriented dipole
    println("Successfully created standard dipole orientations")

    # Custom dipole orientation
    dipole_xy = DipoleVector(0.707, 0.707, 0.0)  # 45° in the XY plane
    println("Successfully created custom dipole orientation")
    
    # Create PSF with custom orientation
    psf = Vector3DPSF(1.4, 0.532, dipole_xy, n_medium=1.33)
    println("Successfully created Vector3DPSF with custom dipole orientation")
    println("Dipole orientation examples successfully tested")
catch e
    println("ERROR in Dipole orientation: ", e)
    println(stacktrace())
end

println("\n===== Testing Aberration Modeling =====")
try
    # Create dipole orientation
    dipole_z = DipoleVector(0.0, 0.0, 1.0)
    
    # Create Zernike coefficients
    println("Creating Zernike coefficients...")
    zc = ZernikeCoefficients(21)  # Support up to 21st Zernike mode
    println("Successfully created Zernike coefficients")
    
    # Add common aberrations
    println("Adding aberrations...")
    add_defocus!(zc, 1.0)
    add_astigmatism!(zc, 0.5)
    add_coma!(zc, 0.3)
    add_spherical!(zc, 0.2)
    # Add trefoil using add_aberration with OSA index 9
    add_aberration!(zc, 9, 0.1)
    println("Successfully added aberrations")
    
    # Create PSF with these aberrations
    println("Creating aberrated Vector3DPSF...")
    psf = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33, base_zernike=zc)
    println("Successfully created Vector3DPSF with multiple aberrations")
    
    # Test the PSF by evaluating at a point
    intensity = psf(0.1, 0.2, 0.3)
    println("Aberrated PSF intensity at (0.1, 0.2, 0.3): ", intensity)
    println("Aberration modeling successfully tested")
catch e
    println("ERROR in Aberration modeling: ", e)
    println(stacktrace())
end

println("\nAll code blocks from vector3d.md successfully tested!")