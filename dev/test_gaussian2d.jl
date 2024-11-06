# dev/test_gaussian2d.jl

using Revise
using MicroscopePSFs
using CairoMakie
using BenchmarkTools
using Printf  # For formatted strings

# Test setup
pixel_size = 0.1  # 100 nm pixels
nx = 20
ny = 12  # Different number of pixels in x and y
σ = 0.13  # 130 nm σ


# Create camera and PSF
x_edges = range(0, nx * pixel_size, nx + 1) |> collect
y_edges = range(0, ny * pixel_size, ny + 1) |> collect
camera = IdealCamera(x_edges, y_edges)
psf = Gaussian2D(σ)

# Centered PSF sampling coordinates (for PSF evaluation)
x = range(-1, 1, length = 100)
y = range(-1, 1, length = 100)

# Compute amplitude values using the built-in 'amplitude' method
amplitude_values = [amplitude(psf, xi, yi) for yi in y, xi in x]

# Break into magnitude and phase
magnitude_values = abs.(amplitude_values)
phase_values = angle.(amplitude_values)

# Compute PSF intensity values
intensity_values = [psf(xi, yi) for yi in y, xi in x]

# Integrated pixels - centered emitter
center_x = (nx - 1) / 2 * pixel_size  # Center of the camera in x
center_y = (ny - 1) / 2 * pixel_size  # Center of the camera in y
emitter_center = Emitter2D(center_x, center_y, 1000.0)

# Use the function 'integrate_pixels_amplitude' from MicroscopePSFs
pixels_center = integrate_pixels(psf, camera, emitter_center)
pixels_center_amp = integrate_pixels_amplitude(psf, camera, emitter_center)

# Integrated pixels - off-center emitter
emitter_off = Emitter2D(center_x + 0.2, center_y - 0.1, 1000.0)
pixels_off = integrate_pixels(psf, camera, emitter_off)
pixels_off_amp = integrate_pixels_amplitude(psf, camera, emitter_off)

# Break integrated amplitude into magnitude and phase
pixels_center_magnitude = abs.(pixels_center_amp)
pixels_center_phase = angle.(pixels_center_amp)
pixels_off_magnitude = abs.(pixels_off_amp)
pixels_off_phase = angle.(pixels_off_amp)

# Create visualization
fig = Figure(size = (1600, 1200))  # Adjusted size for additional plots

# Helper function for consistent axis formatting
function setup_axis!(ax, title)
    ax.title = title
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"
    ax.aspect = DataAspect()
end

# Top row: PSF amplitude magnitude, phase, and intensity (centered at 0,0)
ax1 = Axis(fig[1, 1], width = 300, height = 300)
setup_axis!(ax1, "Amplitude Magnitude (0,0 centered)")
hm1 = heatmap!(ax1, x, y, magnitude_values, colormap = :viridis)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3], width = 300, height = 300)
setup_axis!(ax2, "Amplitude Phase (0,0 centered)")
hm2 = heatmap!(ax2, x, y, phase_values, colormap = :twilight)
Colorbar(fig[1, 4], hm2)

ax_intensity = Axis(fig[1, 5], width = 300, height = 300)
setup_axis!(ax_intensity, "PSF Intensity (0,0 centered)")
hm_intensity = heatmap!(ax_intensity, x, y, intensity_values, colormap = :viridis)
Colorbar(fig[1, 6], hm_intensity)

# Second row: Integrated intensity pixels with emitter coordinates in titles
# Physical coordinates for pixel centers (computed from camera edges)
x_phys = (camera.pixel_edges_x[1:end - 1] + camera.pixel_edges_x[2:end]) / 2
y_phys = (camera.pixel_edges_y[1:end - 1] + camera.pixel_edges_y[2:end]) / 2

ax3 = Axis(fig[2, 1], width = 300, height = 300)
title_center = @sprintf("Centered Emitter at (%.2f, %.2f) μm", emitter_center.x, emitter_center.y)
setup_axis!(ax3, title_center)
ax3.yreversed = true
hm3 = heatmap!(ax3, x_phys, y_phys, pixels_center', colormap = :viridis)
scatter!(ax3, [emitter_center.x], [emitter_center.y],
    color = :red, marker = :cross, markersize = 15)
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3], width = 300, height = 300)
title_off = @sprintf("Off-center Emitter at (%.2f, %.2f) μm", emitter_off.x, emitter_off.y)
setup_axis!(ax4, title_off)
ax4.yreversed = true
hm4 = heatmap!(ax4, x_phys, y_phys, pixels_off', colormap = :viridis)
scatter!(ax4, [emitter_off.x], [emitter_off.y],
    color = :red, marker = :cross, markersize = 15)
Colorbar(fig[2, 4], hm4)

# Third row: Integrated amplitude magnitude and phase
ax5 = Axis(fig[3, 1], width = 300, height = 300)
setup_axis!(ax5, "Integrated Amplitude Magnitude (Centered)")
ax5.yreversed = true
hm5 = heatmap!(ax5, x_phys, y_phys, pixels_center_magnitude', colormap = :viridis)
scatter!(ax5, [emitter_center.x], [emitter_center.y],
    color = :red, marker = :cross, markersize = 15)
Colorbar(fig[3, 2], hm5)

ax6 = Axis(fig[3, 3], width = 300, height = 300)
setup_axis!(ax6, "Integrated Amplitude Phase (Centered)")
ax6.yreversed = true
hm6 = heatmap!(ax6, x_phys, y_phys, pixels_center_phase', colormap = :twilight)
scatter!(ax6, [emitter_center.x], [emitter_center.y],
    color = :red, marker = :cross, markersize = 15)
Colorbar(fig[3, 4], hm6)

# Add coordinate system information
Label(fig[0, :], "Gaussian2D PSF Tests (σ = $σ μm, pixel size = $(pixel_size * 1000) nm)")

# Print validation info
println("\nValidation:")
println("Array dimensions [y, x]: ", size(pixels_center))
println("Expected dimensions [y, x]: [$ny, $nx]")

println("\nCoordinate ranges:")
println("PSF x: ", extrema(x))
println("PSF y: ", extrema(y))
println("Camera x (x_phys): ", extrema(x_phys))
println("Camera y (y_phys): ", extrema(y_phys))

println("\nEmitter positions (physical coordinates):")
println("Centered: (", emitter_center.x, ", ", emitter_center.y, ")")
println("Off-center: (", emitter_off.x, ", ", emitter_off.y, ")")

# Save figure
save("dev/gaussian2d_test.png", fig)

# Return figure handle
fig
