using MicroscopePSFs
using CairoMakie

println("Testing code blocks from interface.md")

# Testing PSF Evaluation code block
println("Testing PSF Evaluation section...")
# 2D evaluation
psf2d = Gaussian2D(0.15)
val_2d = psf2d(0.1, 0.2)
println("2D evaluation: ", val_2d)

# 3D evaluation
psf3d = Scalar3DPSF(1.4, 0.532, 1.518)
val_3d = psf3d(0.1, 0.2, 0.3)
println("3D evaluation: ", val_3d)

# Get complex amplitude
amp_2d = amplitude(psf2d, 0.1, 0.2)
println("2D amplitude: ", amp_2d)
amp_3d = amplitude(psf3d, 0.1, 0.2, 0.3)
println("3D amplitude: ", amp_3d)

# Testing Pixel Integration code block
println("\nTesting Pixel Integration section...")
# Create camera and emitter
# Using constructor with pixel edges directly
pixel_edges_x = collect(0:0.1:2.0)  # Convert to Vector
pixel_edges_y = collect(0:0.1:2.0)  # Convert to Vector
camera = IdealCamera(pixel_edges_x, pixel_edges_y)  # 20x20 pixels, 100nm size
emitter = Emitter2D(1.0, 1.0, 1000.0)               # At (1μm, 1μm) with 1000 photons

# Generate realistic microscope image
pixels = integrate_pixels(psf2d, camera, emitter)
println("Pixel integration dimensions: ", size(pixels))

# For complex amplitude integration
pixels_amplitude = integrate_pixels_amplitude(psf2d, camera, emitter)
println("Amplitude integration dimensions: ", size(pixels_amplitude))

# Testing Creating PSF Instances code block
println("\nTesting Creating PSF Instances section...")
# Create a Gaussian2D PSF with sigma=150nm
psf_gaussian = Gaussian2D(0.15)

# Create an Airy2D PSF with NA=1.4 and wavelength=532nm
psf_airy = Airy2D(1.4, 0.532)

# Create a Vector3D PSF with more parameters
# Create a dipole vector for the orientation (z-axis in this case)
dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Dipole along z-axis
psf_vector = Vector3DPSF(
    1.4,                # Numerical aperture
    0.68,               # Wavelength in microns
    dipole_z,           # Dipole orientation (along optical axis)
    n_medium=1.33,      # Sample medium refractive index
    n_immersion=1.518,  # Immersion medium refractive index
    n_coverslip=1.518   # Coverslip refractive index
)
println("PSF instances created successfully")

# Testing Using the Same Code with Different PSF Models code block
println("\nTesting Using the Same Code with Different PSF Models section...")
function analyze_psf_width(psf::AbstractPSF)
    # Generate intensity profile
    x = range(-1, 1, length=100)
    intensities = [psf(xi, 0.0) for xi in x]
    
    # Calculate FWHM or other properties
    # ...
    
    return "Analysis completed"
end

# Works with any PSF model that implements the 2D interface
results_gaussian = analyze_psf_width(Gaussian2D(0.15))
results_airy = analyze_psf_width(Airy2D(1.4, 0.532))
println("Analysis results: ", results_gaussian, ", ", results_airy)

# Testing Visualizing Different PSF Models code block
println("\nTesting Visualizing Different PSF Models section...")

# Define position and image grid
x = range(-1, 1, length=100)
y = range(-1, 1, length=100)

# Create different PSF models - using only 2D PSFs
psfs = [
    Gaussian2D(0.15),
    Airy2D(1.4, 0.532)
    # Scalar3DPSF only supports 3D interface with (x,y,z) and isn't included here
    # Vector3DPSF only supports 3D interface and also isn't included here
]
titles = ["Gaussian2D", "Airy2D"]

# Compute PSF intensity values
intensity_values = [[psf(xi, yi) for yi in y, xi in x] for psf in psfs]

# Create visualization
fig = Figure(size=(1000, 800))
for (i, img) in enumerate(intensity_values)
    ax = Axis(fig[div(i-1, 2)+1, mod(i-1, 2)+1], 
              aspect=DataAspect(),
              title=titles[i],
              xlabel="x (μm)",
              ylabel="y (μm)")
    hm = heatmap!(ax, x, y, img, colormap=:viridis)
    Colorbar(fig[div(i-1, 2)+1, mod(i-1, 2)+3], hm)
end

# Display or save the figure
#save("interface_psf_visualization.png", fig) # Uncomment to save
println("Visualization created successfully")
display(fig)

println("\nAll code blocks from interface.md tested successfully!")