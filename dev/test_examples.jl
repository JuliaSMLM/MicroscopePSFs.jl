using Pkg 
Pkg.activate("dev")
using Revise
using MicroscopePSFs
using CairoMakie
using Enzyme 

## Render a camera image as a heatmap 

# Create a psf 
psf = AiryPSF(1.4, 0.532)

# Setup camera 
nx, ny, pixel_size = 20, 10, 0.05  # 100 nm pixels
camera = IdealCamera(nx, ny, pixel_size)

# Render an image with a single emitter
x_microns, y_microns, photons = 0.5, 0.25, 1000.0  # Emitter position and brightness
image = integrate_pixels(psf, camera, Emitter2D(x_microns, y_microns, photons))

# Show image as a heatmap with pixel (1,1) at top left and y axis descending
fig = Figure()
ax = Axis(fig[1, 1], yreversed=true, aspect=DataAspect())
heatmap!(ax, image', colormap=:inferno)
display(fig)
save("camera_heatmap.png", fig)


## Generate a simulation with multiple emitters

# Create a psf with astigmatism
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5  # Add vertical astigmatism 
psf = ScalarPSF(1.4, 0.532, 1.52; zernike_coeffs=zc)

# Create SplinePSF for speed. This is the slowest part of the simulation.
xy_sampling, z_sampling = 0.05, 0.1
x_range = y_range = -1.0:xy_sampling:1.0
z_range = -1.0:xy_sampling:1.0
psf_spline = SplinePSF(psf, x_range, y_range, z_range)

# Setup camera
nx, ny, pixel_size = 128, 128, 0.1  # 100 nm pixels
camera = IdealCamera(nx, ny, pixel_size)

# Create emitters
n_emitters = 5.0 * (nx * ny * pixel_size^2)  # 5 μm^-2 of emitter density
emitters = [Emitter3D(nx * pixel_size * rand(), ny * pixel_size * rand(), rand()-0.5, 1000.0) 
    for _ in 1:n_emitters]

# Render the image using region of support option
image = integrate_pixels(psf_spline, camera, emitters; support=0.5)
  

fig = Figure()
ax = Axis(fig[1, 1], yreversed=true, aspect=DataAspect())
heatmap!(ax, image', colormap=:inferno)
hidedecorations!(ax)
display(fig)
save("camera_multi_emitters.png", fig)



## Model a PSF measurement with stage movement (as opposed to emitter movement)

# Create a Vector PSF with x-dipole 
dipole = DipoleVector(1.0, 0.0, 0.0)  # x-dipole
emitter = Emitter3D(0.5, 0.5, 0.0, 1000.0)  # Emitter on coverslip

# Setup camera
nx, ny, pixel_size = 20, 20, 0.05  # 100 nm pixels
camera = IdealCamera(nx, ny, pixel_size)

z_stack = Vector{Array{Float64}}()

# Loop over z_stage positions
stage_positions = -1.0:0.5:1.0  # Z stage positions in microns
for z_stage in stage_positions
    # Create a VectorPSF with the current z_stage position
    psf = VectorPSF(1.4, 0.690, dipole; z_stage=z_stage)
    
    # Integrate over the camera pixels
    image = integrate_pixels(psf, camera, emitter)
    
    # Store the image in the z_stack
    push!(z_stack, image)
end

fig = Figure()
for (i, image) in enumerate(z_stack)
    ax = Axis(fig[1, i], yreversed=true, aspect=DataAspect(), title="$(stage_positions[i]) μm")
    heatmap!(ax, image', colormap=:inferno)
    hidedecorations!(ax)
end
display(fig)
save("camera_z_stack_stage.png", fig)

## Example of using Enzyme to differentiate a PSF 

# Camera setup
nx, ny, pixel_size = 20, 20, 0.1  # 100 nm pixels
camera = IdealCamera(nx, ny, pixel_size)
camera = IdealCamera(nx, ny, pixel_size)

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5  # Add vertical astigmatism
psf = ScalarPSF(1.2, 0.6, 1.33; zernike_coeffs=zc)

# Calc derivative image 
x_range = y_range = range(-0.5, 0.5, 40)  # PSF field coordinates

# Calculate the PSF and its derivatives using Enzyme
dx_values = [Enzyme.autodiff(Enzyme.Reverse, x -> psf(x, y, 0.0), Active, Active(x))[1][1] for x in x_range, y in y_range]
dy_values = [Enzyme.autodiff(Enzyme.Reverse, y -> psf(x, y, 0.0), Active, Active(y))[1][1] for x in x_range, y in y_range]
dz_values = [Enzyme.autodiff(Enzyme.Reverse, z -> psf(x, y, z), Active, Active(0.0))[1][1] for x in x_range, y in y_range]
psf_arr = [psf(x, y, 0.0) for x in x_range, y in y_range]

# Plot with CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], title = "PSF", yreversed=true, aspect=DataAspect())
heatmap!(ax, psf_arr')
ax2 = Axis(fig[1, 2], title = "dPSF/dx", yreversed=true, aspect=DataAspect())
heatmap!(ax2, dx_values')
ax3 = Axis(fig[1, 3], title = "dPSF/dy", yreversed=true, aspect=DataAspect())
heatmap!(ax3, dy_values')
ax4 = Axis(fig[1, 4], title = "dPSF/dz", yreversed=true, aspect=DataAspect())
heatmap!(ax4, dz_values')
display(fig)
save("psf_derivatives.png", fig)

# Now try with integrating the PSF. 

# Setup camera
nx, ny, pixel_size = 20, 20, 0.1  # 100 nm pixels
camera = IdealCamera(nx, ny, pixel_size)

x_emitter = y_emitter = nx / 2 * pixel_size
z_emitter, photons = 0.0, 1000.0

# Use jacobian since our output is a 2D image
dx_image = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
    x -> integrate_pixels(psf, camera, Emitter3D(x, y_emitter, z_emitter, photons)), 
    x_emitter)[1]
dy_image = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
    y -> integrate_pixels(psf, camera, Emitter3D(x_emitter, y, z_emitter, photons)), 
    y_emitter)[1]
dz_image = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
    z -> integrate_pixels(psf, camera, Emitter3D(x_emitter, y_emitter, z, photons)), 
    z_emitter)[1]
psf_image = integrate_pixels(psf, camera, Emitter3D(x_emitter, y_emitter, z_emitter, photons))

# Plot with CairoMakie
fig2 = Figure()
ax5 = Axis(fig2[1, 1], title = "Integrated PSF", yreversed=true, aspect=DataAspect())
heatmap!(ax5, psf_image')
ax6 = Axis(fig2[1, 2], title = "dPSF/dx", yreversed=true, aspect=DataAspect())
heatmap!(ax6, dx_image')
ax7 = Axis(fig2[1, 3], title = "dPSF/dy", yreversed=true, aspect=DataAspect())
heatmap!(ax7, dy_image')
ax8 = Axis(fig2[1, 4], title = "dPSF/dz", yreversed=true, aspect=DataAspect())
heatmap!(ax8, dz_image')
display(fig2)
save("psf_integrated_derivatives.png", fig2)





