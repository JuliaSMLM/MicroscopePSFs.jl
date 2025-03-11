# dev/test_vector3d.jl
using Pkg 
Pkg.activate("dev")

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
using SMLMData

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n_medium = 1.33   # Sample medium refractive index (water)
n_coverslip = 1.52  # Coverslip refractive index (glass)
n_immersion = 1.52  # Immersion medium refractive index (oil)
focal_z = 1.0  # Nominal focal plane position

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(0:nx, 0:ny, pixel_size)
cam_size = nx * pixel_size

# Create different dipole orientations to test
dipoles = [
    DipoleVector(1.0, 0.0, 0.0),    # x-oriented
    DipoleVector(0.0, 1.0, 0.0),    # y-oriented
    DipoleVector(0.0, 0.0, 1.0),    # z-oriented
    DipoleVector(1.0, 1.0, 1.0),    # xyz-oriented (will be normalized)
]

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.0  # Add vertical astigmatism for visualization

# Create Vector3DPSF with the first dipole orientation (x-oriented)
psf = Vector3DPSF(
    na, λ, dipoles[1]; 
    base_zernike=zc,
    n_medium=n_medium, 
    n_coverslip=n_coverslip, 
    n_immersion=n_immersion,
    focal_z=focal_z
)

# Check normalization over camera
emitter_center = Emitter3D(nx/2 * pixel_size, ny/2 * pixel_size, 0.0, 0.0)
psf_normalized = integrate_pixels(psf, camera, emitter_center)
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

# Get camera coordinates for display
x_cam = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) ./ 2
y_cam = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) ./ 2

# Row 1: Pupil field components - Ex and Ey amplitude
ax_ex_amp = Axis(fig[1, 1], width=225, height=225)
setup_axis!(ax_ex_amp, "Ex Pupil Amplitude")
hm_ex_amp = heatmap!(ax_ex_amp, x, y, abs.(psf.vector_pupils.Ex.field), colormap=:viridis)
Colorbar(fig[1, 2], hm_ex_amp)

ax_ey_amp = Axis(fig[1, 3], width=225, height=225)
setup_axis!(ax_ey_amp, "Ey Pupil Amplitude")
hm_ey_amp = heatmap!(ax_ey_amp, x, y, abs.(psf.vector_pupils.Ey.field), colormap=:viridis)
Colorbar(fig[1, 4], hm_ey_amp)

# Row 2: Different dipole orientations (at z=0)
for (i, dipole) in enumerate(dipoles)
    ax = Axis(fig[2, i], width=225, height=225)
    setup_axis!(ax, "Dipole: $(round.(typeof(dipole.px).([dipole.px, dipole.py, dipole.pz]), digits=1))")
    
    # Create PSF with this dipole orientation
    @time dipole_psf = Vector3DPSF(
        na, λ, dipole; 
        base_zernike=zc,
        n_medium=n_medium, 
        n_coverslip=n_coverslip, 
        n_immersion=n_immersion,
        focal_z=focal_z
    )
    
    # Compute intensity
    intensity = [dipole_psf(xi, yi, 0.0) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, intensity, colormap=:viridis)
    
    if i == length(dipoles)
        Colorbar(fig[2, i+1], hm)
    end
end

# Row 3: PSF field through focus (with the primary dipole orientation)
for (i, z) in enumerate(z_planes)
    ax = Axis(fig[3, i], width=225, height=225)
    setup_axis!(ax, @sprintf("z = %.1f μm", z))
    intensity = [psf(xi, yi, z) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, intensity, colormap=:viridis)
    if i == length(z_planes)
        Colorbar(fig[3, i+1], hm)
    end
end

# Row 4: Integrated pixels for different emitter positions
for (i, emitter) in enumerate(emitters)
    ax = Axis(fig[4, i], width=225, height=225)
    setup_axis!(ax, @sprintf("(%.1f, %.1f, %.1f)μm", emitter.x, emitter.y, emitter.z))
    ax.yreversed = true
    
    pixels = integrate_pixels(psf, camera, emitter)
    hm = heatmap!(ax, x_cam, y_cam, pixels', colormap=:viridis)
    scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
    
    if i == length(emitters)
        Colorbar(fig[4, i+1], hm)
    end
end

# Row 5: Ex and Ey field components (for one emitter position)
emitter_focus = emitters[1]  # Use the in-focus emitter
field_components = integrate_pixels_amplitude(psf, camera, emitter_focus)

ax_ex = Axis(fig[5, 1], width=225, height=225)
setup_axis!(ax_ex, "Ex amplitude")
ax_ex.yreversed = true
hm_ex = heatmap!(ax_ex, x_cam, y_cam, abs.(field_components[:,:,1])', colormap=:viridis)
Colorbar(fig[5, 2], hm_ex)

ax_ey = Axis(fig[5, 3], width=225, height=225)
setup_axis!(ax_ey, "Ey amplitude")
ax_ey.yreversed = true
hm_ey = heatmap!(ax_ey, x_cam, y_cam, abs.(field_components[:,:,2])', colormap=:viridis)
Colorbar(fig[5, 4], hm_ey)



# Add labels and title
Label(fig[1, 0], "Pupil\nComponents", rotation=π/2)
Label(fig[2, 0], "Dipole\nOrientation", rotation=π/2)
Label(fig[3, 0], "PSF\nFocus", rotation=π/2)
Label(fig[4, 0], "Camera\nPixels", rotation=π/2)

title_text = @sprintf("Vector3D PSF with Astigmatism\n(NA=%.1f, λ=%.3f μm, n₁=%.2f, n₂=%.2f, n₃=%.2f)", 
                      na, λ, n_medium, n_coverslip, n_immersion)
Label(fig[0, :], title_text, fontsize=20, padding=(0,0,20,0))

# Adjust layout
rowgap!(fig.layout, 40)
colgap!(fig.layout, 40)

# Save figure
save("dev/vector3d_test.png", fig)

# Print parameters
println("\nPSF parameters:")
println("NA = $na, λ = $λ μm")
println("Medium index = $n_medium")
println("Coverslip index = $n_coverslip")
println("Immersion index = $n_immersion")
println("Focal plane = $focal_z μm")
println("Astigmatism: 1 wave at 0 degrees")
println("Field size: $(cam_size) × $(cam_size) μm")
println("Pixel size: $(pixel_size) μm")
println("Z planes: ", z_planes, " μm")
println("Dipole orientations:")
for (i, dipole) in enumerate(dipoles)
    println("  $i: ($(dipole.px), $(dipole.py), $(dipole.pz))")
end

fig