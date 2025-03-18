# multi_emitter_demo.jl

using Pkg
Pkg.activate("dev")

using Revise
using MicroscopePSFs
using CairoMakie


# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 60           # Width
ny = 40           # Height
camera = IdealCamera(nx, ny, pixel_size)

# Create a PSF
psf = Scalar3DPSF(na, λ, n)

# Create multiple emitters in a grid pattern with slight variations
emitters = AbstractEmitter[]
photon_range = 800.0:200.0:1800.0  # Different brightness levels
center_x, center_y = nx / 2 * pixel_size, ny / 2 * pixel_size  # Center of the camera

# Add emitters in a grid pattern with different z positions
for i in 1:3
    for j in 1:3
        x = center_x + (i - 2) * 0.3  # 300 nm spacing
        y = center_y + (j - 2) * 0.3
        z = 0.2 * ((i - 2) + (j - 2))   # Vary z position for visual effect
        photons = photon_range[i+j-1]  # Vary brightness
        push!(emitters, Emitter3D(x, y, z, photons))
    end
end

# Add a special emitter offset from the grid
push!(emitters, Emitter3D(center_x + 0.8, center_y - 0.5, -0.3, 2000.0))

# Single emitter for comparison
single_emitter = emitters[5]  # Middle emitter

# Integrate pixels for all emitters at once
println("Integrating $(length(emitters)) emitters with finite support region...")
@time multi_image = integrate_pixels(psf, camera, emitters; support=0.5)

# Integrate the single middle emitter for comparison
println("Integrating single emitter...")
@time single_image = integrate_pixels(psf, camera, single_emitter)


# Create visualization
fig = Figure(size=(900, 400))

# Helper function for axis setup
function setup_axis!(ax, title)
    ax.title = title
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"
    ax.aspect = DataAspect()
    ax.yreversed = true  # To match camera coordinates (0,0 is top-left)
end

# Get camera coordinates for display
x_cam = (camera.pixel_edges_x[1:end-1] + camera.pixel_edges_x[2:end]) ./ 2
y_cam = (camera.pixel_edges_y[1:end-1] + camera.pixel_edges_y[2:end]) ./ 2

# Plot multi-emitter image
ax1 = Axis(fig[1, 1])
setup_axis!(ax1, "Multiple Emitters ($(length(emitters)))")
hm1 = heatmap!(ax1, x_cam, y_cam, multi_image', colormap=:viridis)
# Mark emitter positions
for emitter in emitters
    scatter!(ax1, [emitter.x], [emitter.y], color=:red, marker=:circle,
        markersize=10 * (emitter.photons / 2000.0)^0.5)  # Size by brightness
end
Colorbar(fig[1, 2], hm1)

# Plot single emitter image
ax2 = Axis(fig[1, 3])
setup_axis!(ax2, "Single Emitter")
hm2 = heatmap!(ax2, x_cam, y_cam, single_image', colormap=:viridis)
scatter!(ax2, [single_emitter.x], [single_emitter.y], color=:red, marker=:circle, markersize=10)
Colorbar(fig[1, 4], hm2)

# Save figure
save("multi_emitter_demo.png", fig)

println("Total photons in multi-emitter image: $(sum(multi_image))")
println("Total photons in single emitter image: $(sum(single_image))")
println("Ratio of sum to emitted photons: $(sum(multi_image) / sum(e.photons for e in emitters))")

# Display the figure
fig