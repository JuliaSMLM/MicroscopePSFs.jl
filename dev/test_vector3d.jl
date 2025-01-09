# dev/test_vector3d.jl

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
using SMLMData

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index of immersion medium
n_med = 1.33  # Refractive index of medium (water)

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(1:nx, 1:ny, pixel_size)
cam_size = nx * pixel_size

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
# Add astigmatism
# zc.phase[6] = 1.0

# Create base pupil function
base_pupil = PupilFunction(na, λ, n, zc)

# Create Vector PSF
psf = Vector3DPupilPSF(na, λ, base_pupil; n_medium=n_med)

# Evaluation coordinates
x = y = range(-1, 1, 100)  # PSF field coordinates
z_planes = [-1.0, -0.5, 0.0, 0.5, 1.0]  # z positions in microns

# Different dipole orientations to test
dipoles = [
    DipoleVector(1.0, 0.0, 0.0),  # x-oriented
    DipoleVector(0.0, 1.0, 0.0),  # y-oriented
    DipoleVector(0.0, 0.0, 1.0)  # z-oriented
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

# Row 1: Pupil amplitude and phase for x-dipole
vpupil = psf.pupil
Ex_field = vpupil.pupils[1].field  # Ex from x-dipole
Ey_field = vpupil.pupils[4].field  # Ey from x-dipole

ax_ex = Axis(fig[1, 1], width=225, height=225)
setup_axis!(ax_ex, "Ex Amplitude (x-dipole)")
hm_ex = heatmap!(ax_ex, x, y, abs.(Ex_field), colormap=:viridis)
Colorbar(fig[1, 2], hm_ex)

ax_ey = Axis(fig[1, 3], width=225, height=225)
setup_axis!(ax_ey, "Ey Amplitude (x-dipole)")
hm_ey = heatmap!(ax_ey, x, y, abs.(Ey_field), colormap=:viridis)
Colorbar(fig[1, 4], hm_ey)

# Row 2: PSF through focus for x-dipole
for (i, z) in enumerate(z_planes)
    ax = Axis(fig[2, i], width=225, height=225)
    setup_axis!(ax, @sprintf("z = %.1f μm", z))
    intensity = [psf(xi, yi, z, dipoles[1]) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, intensity, colormap=:viridis)
    if i == length(z_planes)
        Colorbar(fig[2, i+1], hm)
    end
end

# Row 3: PSF for different dipole orientations at z=0
for (i, dipole) in enumerate(dipoles)
    ax = Axis(fig[3, i], width=225, height=225)
    setup_axis!(ax, @sprintf("(%.1f, %.1f, %.1f)", dipole.px, dipole.py, dipole.pz))
    intensity = [psf(xi, yi, 0.0, dipole) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, intensity, colormap=:viridis)
    if i == length(dipoles)
        Colorbar(fig[3, i+1], hm)
    end
end

# Row 4: Camera pixels for isotropic emitters at different z
z_test = [-0.5, 0.0, 0.5]
for (i, z) in enumerate(z_test)
    # Get camera coordinates for display
    x_cam = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) ./ 2
    y_cam = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) ./ 2
    
    ax = Axis(fig[4, i], width=225, height=225)
    setup_axis!(ax, @sprintf("z = %.1f μm", z))
    ax.yreversed = true
    
    emitter = Emitter3D(1.0, 0.5, z, 1000.0)
    pixels = integrate_pixels(psf, camera, emitter)
    hm = heatmap!(ax, x_cam, y_cam, pixels', colormap=:viridis)
    scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    
    if i == length(z_test)
        Colorbar(fig[4, i+1], hm)
    end
end

# Add labels
Label(fig[1, 0], "Pupil\nFields", rotation=π/2)
Label(fig[2, 0], "X-dipole\nPSF", rotation=π/2)
Label(fig[3, 0], "Dipole\nOrientations", rotation=π/2)
Label(fig[4, 0], "Camera\nPixels", rotation=π/2)

Label(fig[0, :], "Vector3D PSF with Astigmatism\n(NA=$na, λ=$λ μm, n=$n, n_medium=$n_med)", 
      fontsize=20, padding=(0,0,20,0))

# Adjust layout
rowgap!(fig.layout, 40)
colgap!(fig.layout, 40)

# Save figure
save("dev/vector3d_test.png", fig)

# Print parameters
println("\nPSF parameters:")
println("NA = $na, λ = $λ μm")
println("n_immersion = $n, n_medium = $n_med")
println("Astigmatism: 1 wave")
println("Field size: $(cam_size) × $(cam_size) μm")
println("Pixel size: $(pixel_size) μm")
println("Z planes: ", z_planes, " μm")

fig