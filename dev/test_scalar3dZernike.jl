# dev/test_scalar3d.jl

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
using SMLMData

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(1:nx, 1:ny, pixel_size)
cam_size = nx * pixel_size

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
# add_astigmatism!(zc, 2.0, 0.0)  # 2 waves of astigmatism at 0 degrees
psf = Scalar3DPSF(na, λ, n, coeffs=zc)
zc.phase[4] = 1.0  # Add phase shift for visualization

# Check normalization over camera
psf_normalized = integrate_pixels(psf, camera, Emitter3D(nx/2 * pixel_size, ny/2 * pixel_size, 0.0, 1.0))
println("PSF normalization: ", sum(psf_normalized))

# Evaluation coordinates
x = y = range(-1, 1, 100)  # PSF field coordinates
z_planes = [-1.0, -0.5, 0.0, 0.5, 1.0]  # z positions in microns

# Emitter positions to test
emitters = [
    Emitter3D(1.0, 0.5, 0.0, 1000.0),    # center, in focus
    Emitter3D(1.5, 0.3, 0.5, 1000.0),    # right side, above focus
    Emitter3D(0.5, 0.7, -0.5, 1000.0),   # left side, below focus
    Emitter3D(1.2, 0.4, 1.0, 1000.0),    # off center, far from focus
]


# Create visualization
fig = Figure(
    size = (2000, 1600),
    backgroundcolor = :white,
    figure_padding = (40, 40, 40, 40)
)

# Helper function for axis setup
function setup_axis!(ax, title)
    ax.title = title
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"
    ax.aspect = DataAspect()
end

# Get pupil field for display
pupil = PupilFunction(na, λ, n, zc, grid_size=100)

# Row 1: Pupil amplitude and phase
ax_amp = Axis(fig[1, 1], width=225, height=225)
setup_axis!(ax_amp, "Pupil Amplitude")
hm_amp = heatmap!(ax_amp, x, y, abs.(pupil.field), colormap=:viridis)
Colorbar(fig[1, 2], hm_amp)

ax_phase = Axis(fig[1, 3], width=225, height=225)
setup_axis!(ax_phase, "Pupil Phase")
hm_phase = heatmap!(ax_phase, x, y, angle.(pupil.field), colormap=:twilight)
Colorbar(fig[1, 4], hm_phase)

# Row 2: PSF field through focus
for (i, z) in enumerate(z_planes)
    ax = Axis(fig[2, i], width=225, height=225)
    setup_axis!(ax, @sprintf("z = %.1f μm", z))
    intensity = [psf(xi, yi, z) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, intensity, colormap=:viridis)
    if i == length(z_planes)
        Colorbar(fig[2, i+1], hm)
    end
end

# Rows 3-4: Integrated pixels for different emitter positions
for (i, emitter) in enumerate(emitters)
    # Get camera coordinates for display
    x_cam = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) ./ 2
    y_cam = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) ./ 2
    
    ax = Axis(fig[3, i], width=225, height=225)
    setup_axis!(ax, @sprintf("(%.1f, %.1f, %.1f)μm", emitter.x, emitter.y, emitter.z))
    ax.yreversed = true
    
    pixels = integrate_pixels(psf, camera, emitter)
    hm = heatmap!(ax, x_cam, y_cam, pixels', colormap=:viridis)
    scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    
    if i == length(emitters)
        Colorbar(fig[3, i+1], hm)
    end
end

# Add labels and title
Label(fig[1, 0], "Pupil", rotation=π/2)
Label(fig[2, 0], "PSF\nField", rotation=π/2)
Label(fig[3, 0], "Camera\nPixels", rotation=π/2)

Label(fig[0, :], "Scalar3D PSF with Astigmatism\n(NA=$na, λ=$λ μm, n=$n)", 
      fontsize=20, padding=(0,0,20,0))

# Adjust layout
rowgap!(fig.layout, 40)
colgap!(fig.layout, 40)

# Save figure
save("dev/scalar3d_test.png", fig)

# Print parameters
println("\nPSF parameters:")
println("NA = $na, λ = $λ μm, n = $n")
println("Astigmatism: 2 waves at 0 degrees")
println("Field size: $(cam_size) × $(cam_size) μm")
println("Pixel size: $(pixel_size) μm")
println("Z planes: ", z_planes, " μm")

fig