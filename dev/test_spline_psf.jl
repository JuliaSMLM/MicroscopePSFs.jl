# dev/test_spline_psf.jl

using Pkg
Pkg.activate("dev")
using Revise
using MicroscopePSFs
using CairoMakie
using BenchmarkTools
using Printf
using Statistics

# Test setup
pixel_size = 0.1  # 100 nm pixels
nx = 20
ny = 20
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index

# Create camera and original PSF
x_edges = range(0, nx * pixel_size, nx + 1) |> collect
y_edges = range(0, ny * pixel_size, ny + 1) |> collect
camera = IdealCamera(x_edges, y_edges)

# Create a 3D PSF to be approximated
original_psf = Scalar3DPSF(na, λ, n)

# Create SplinePSF from the original PSF
z_range = range(-2.0, 2.0, length=21)  # 21 z-planes, 4 µm range
sampling = 2  # Use 2× oversampling

# Time the SplinePSF creation
@time spline_psf = SplinePSF(
    original_psf,
    x_range = range(-1.5, 1.5, length=31),  # 31 x samples over 3 µm
    y_range = range(-1.5, 1.5, length=31),  # 31 y samples over 3 µm
    z_range = z_range,                      # 21 z samples over 4 µm
    sampling = sampling                      # 2× oversampling
)

# Print PSF information
println("\nCreated SplinePSF from Scalar3DPSF:")
println("  Original PSF: Scalar3DPSF(na=$na, λ=$λ, n=$n)")
println("  Spline grid dimensions: ", size(spline_psf.coefficients)[1:3])
println("  Spline grid x range: ", extrema(spline_psf.x_knots), " µm")
println("  Spline grid y range: ", extrema(spline_psf.y_knots), " µm")
println("  Spline grid z range: ", extrema(spline_psf.z_knots), " µm")

# Create a fine evaluation grid for visualization
x = y = range(-1.0, 1.0, length=101)

# Select a few z planes to visualize
z_planes = [-1.0, -0.5, 0.0, 0.5, 1.0]

# Calculate timing for evaluation
function evaluate_all_points(psf, x, y, z)
    result = zeros(length(y), length(x), length(z))
    for (iz, z_val) in enumerate(z)
        for (iy, y_val) in enumerate(y)
            for (ix, x_val) in enumerate(x)
                result[iy, ix, iz] = psf(x_val, y_val, z_val)
            end
        end
    end
    return result
end

# Time each PSF evaluation (coarser grid for faster timing)
x_coarse = y_coarse = range(-1.0, 1.0, length=41)
z_coarse = range(-1.0, 1.0, length=5)

println("\nTimings:")
print("Original PSF evaluation: ")
@time original_values = evaluate_all_points(original_psf, x_coarse, y_coarse, z_coarse)
print("Spline PSF evaluation: ")
@time spline_values = evaluate_all_points(spline_psf, x_coarse, y_coarse, z_coarse)

# Calculate relative error
rel_error = abs.(original_values - spline_values) ./ max.(1e-10, abs.(original_values))
mean_rel_error = mean(rel_error)
max_rel_error = maximum(rel_error)
println("Mean relative error: $(mean_rel_error)")
println("Max relative error: $(max_rel_error)")

# Function to compute PSF values on grid
function compute_psf_grid(psf, x, y, z)
    return [psf(x_val, y_val, z_val) for y_val in y, x_val in x, z_val in z]
end

# Compute PSF values for visualization
original_psf_values = compute_psf_grid(original_psf, x, y, z_planes)
spline_psf_values = compute_psf_grid(spline_psf, x, y, z_planes)

# Compute absolute error
abs_error = abs.(original_psf_values - spline_psf_values)

# Create a visualization
fig = Figure(size=(1200, 1200))

# Helper function for consistent axis formatting
function setup_axis!(ax, title)
    ax.title = title
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"
    ax.aspect = DataAspect()
end

# Color range for intensity plots
intensity_colorrange = (0, maximum(original_psf_values) * 0.8)
error_colorrange = (0, max(1e-6, maximum(abs_error) * 0.8))

# Plot PSFs at different z planes
for (i, z) in enumerate(z_planes)
    # Original PSF
    ax1 = Axis(fig[i, 1])
    setup_axis!(ax1, @sprintf("Original PSF, z = %.1f μm", z))
    hm1 = heatmap!(ax1, x, y, original_psf_values[:, :, i],
                  colormap=:viridis, colorrange=intensity_colorrange)
    
    # Spline PSF
    ax2 = Axis(fig[i, 2])
    setup_axis!(ax2, @sprintf("Spline PSF, z = %.1f μm", z))
    hm2 = heatmap!(ax2, x, y, spline_psf_values[:, :, i],
                  colormap=:viridis, colorrange=intensity_colorrange)
    
    # Error
    ax3 = Axis(fig[i, 3])
    setup_axis!(ax3, @sprintf("Absolute Error, z = %.1f μm", z))
    hm3 = heatmap!(ax3, x, y, abs_error[:, :, i],
                  colormap=:viridis, colorrange=error_colorrange)
    
    # Colorbar for the row
    if i == 1
        Colorbar(fig[1, 4], hm2, label="Intensity")
        Colorbar(fig[1, 5], hm3, label="Error")
    end
end

# Add a row for cross-sections
z_idx = 3  # Use z=0 plane
x_idx = 51  # Middle of the x range

ax_cross_x = Axis(fig[6, 1:2])
ax_cross_x.xlabel = "x (μm)"
ax_cross_x.ylabel = "Intensity"
ax_cross_x.title = "Cross-section at y=0, z=0"
lines!(ax_cross_x, x, original_psf_values[x_idx, :, z_idx], label="Original", linewidth=2)
lines!(ax_cross_x, x, spline_psf_values[x_idx, :, z_idx], label="Spline", linewidth=2, linestyle=:dash)
axislegend(ax_cross_x)

ax_cross_error = Axis(fig[6, 3:5])
ax_cross_error.xlabel = "x (μm)"
ax_cross_error.ylabel = "Error"
ax_cross_error.title = "Absolute Error at y=0, z=0"
lines!(ax_cross_error, x, abs_error[x_idx, :, z_idx], linewidth=2, color=:red)

# Now test with camera integration
println("\nCamera Integration Test:")

# Create emitters at different z positions
center_x = (nx - 1) / 2 * pixel_size
center_y = (ny - 1) / 2 * pixel_size
emitter_z0 = Emitter3D(center_x, center_y, 0.0, 1000.0)
emitter_z1 = Emitter3D(center_x, center_y, 1.0, 1000.0)

# Time camera integration
print("Original PSF camera integration: ")
@time pixels_orig_z0 = integrate_pixels(original_psf, camera, emitter_z0, sampling=2)
@time pixels_orig_z1 = integrate_pixels(original_psf, camera, emitter_z1, sampling=2)

print("Spline PSF camera integration: ")
@time pixels_spline_z0 = integrate_pixels(spline_psf, camera, emitter_z0, sampling=2)
@time pixels_spline_z1 = integrate_pixels(spline_psf, camera, emitter_z1, sampling=2)

# Calculate errors for camera integration
camera_error_z0 = abs.(pixels_orig_z0 - pixels_spline_z0)
camera_error_z1 = abs.(pixels_orig_z1 - pixels_spline_z1)
max_camera_error = max(maximum(camera_error_z0), maximum(camera_error_z1))
println("Maximum camera integration error: $max_camera_error")

# Add rows for camera integration comparison
x_cam = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) ./ 2
y_cam = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) ./ 2

# Row for z=0 emitter
ax_cam_orig_z0 = Axis(fig[7, 1])
setup_axis!(ax_cam_orig_z0, "Original PSF, z = 0.0 μm")
ax_cam_orig_z0.yreversed = true
hm_cam1 = heatmap!(ax_cam_orig_z0, x_cam, y_cam, pixels_orig_z0',
                  colormap=:viridis, colorrange=(0, 0.05))

ax_cam_spline_z0 = Axis(fig[7, 2])
setup_axis!(ax_cam_spline_z0, "Spline PSF, z = 0.0 μm")
ax_cam_spline_z0.yreversed = true
heatmap!(ax_cam_spline_z0, x_cam, y_cam, pixels_spline_z0',
        colormap=:viridis, colorrange=(0, 0.05))

ax_cam_error_z0 = Axis(fig[7, 3])
setup_axis!(ax_cam_error_z0, "Integration Error, z = 0.0 μm")
ax_cam_error_z0.yreversed = true
hm_cam_err = heatmap!(ax_cam_error_z0, x_cam, y_cam, camera_error_z0',
                     colormap=:viridis, colorrange=(0, max_camera_error))

# Row for z=1.0 emitter
ax_cam_orig_z1 = Axis(fig[8, 1])
setup_axis!(ax_cam_orig_z1, "Original PSF, z = 1.0 μm")
ax_cam_orig_z1.yreversed = true
heatmap!(ax_cam_orig_z1, x_cam, y_cam, pixels_orig_z1',
        colormap=:viridis, colorrange=(0, 0.05))

ax_cam_spline_z1 = Axis(fig[8, 2])
setup_axis!(ax_cam_spline_z1, "Spline PSF, z = 1.0 μm")
ax_cam_spline_z1.yreversed = true
heatmap!(ax_cam_spline_z1, x_cam, y_cam, pixels_spline_z1',
        colormap=:viridis, colorrange=(0, 0.05))

ax_cam_error_z1 = Axis(fig[8, 3])
setup_axis!(ax_cam_error_z1, "Integration Error, z = 1.0 μm")
ax_cam_error_z1.yreversed = true
heatmap!(ax_cam_error_z1, x_cam, y_cam, camera_error_z1',
        colormap=:viridis, colorrange=(0, max_camera_error))

# Add colorbars for camera integration
Colorbar(fig[7, 4], hm_cam1, label="Intensity")
Colorbar(fig[7, 5], hm_cam_err, label="Error")

# Add title
Label(fig[0, 1:5], text="SplinePSF Test - Comparison with Scalar3DPSF (sampling=$sampling×)",
     fontsize=20)

# Save figure
save("dev/spline_psf_test.png", fig)

println("\nTest complete. Figure saved to dev/spline_psf_test.png")

# Return figure for display
fig