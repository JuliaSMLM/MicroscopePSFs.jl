using Pkg
Pkg.activate("dev")
using Revise
using MicroscopePSFs
using CairoMakie
using LinearAlgebra
using Statistics
using BenchmarkTools

#= PART 1: BASIC USAGE OF SPLINEPSF =#

# Create a ScalarPSF
scalar_psf = ScalarPSF(1.4, 0.532, 1.518)

# Sample to create a SplinePSF
x_range = y_range = range(-2.0, 2.0, step=0.05)  # 41×41 lateral grid, 4μm range
z_range = range(-1.0, 1.0, step=0.1)            # 21 z-planes, 2μm range

println("Type of x_range: ", typeof(x_range))
println("Type of y_range: ", typeof(y_range))
println("Type of z_range: ", typeof(z_range))

# Update the SplinePSF creation to use the positional constructor
println("Creating SplinePSF from ScalarPSF...")
spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)

println("Original PSF type: $(typeof(scalar_psf))")
println("Spline PSF type: $(typeof(spline_psf))")

# Check normalization over camera
# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(1:nx, 1:ny, pixel_size)
cam_size = nx * pixel_size

# Check normalization over camera
pixels_spline = integrate_pixels(spline_psf, camera, Emitter3D(nx / 2 * pixel_size, ny / 2 * pixel_size, 0.0, 1.0))
pixels_scalar = integrate_pixels(scalar_psf, camera, Emitter3D(nx / 2 * pixel_size, ny / 2 * pixel_size, 0.0, 1.0))
println("SplinePSF normalization: ", sum(pixels_spline))
println("ScalarPSF normalization: ", sum(pixels_scalar))


# Create grids for visualization
xgrid = ygrid = range(-1.0, 1.0, length=101)  # Finer grid for visualization
zgrid = [0.0, 0.2, 0.4, 0.6]  # Few z-planes to visualize

# Compare PSFs at different z-planes
println("Comparing PSFs at different z-planes...")

# Function to find global intensity range
function find_intensity_range(scalar_psf, spline_psf, xgrid, ygrid, zgrid)
    local_min = Inf
    local_max = -Inf

    for z in zgrid
        original = [scalar_psf(x, y, z) for y in ygrid, x in xgrid]
        spline = [spline_psf(x, y, z) for y in ygrid, x in xgrid]

        local_min = min(local_min, minimum(original), minimum(spline))
        local_max = max(local_max, maximum(original), maximum(spline))
    end

    return local_min, local_max
end

# Function to find maximum difference
function find_max_difference(scalar_psf, spline_psf, xgrid, ygrid, zgrid)
    max_diff = 0.0

    for z in zgrid
        original = [scalar_psf(x, y, z) for y in ygrid, x in xgrid]
        spline = [spline_psf(x, y, z) for y in ygrid, x in xgrid]
        diff = original - spline
        max_diff = max(max_diff, maximum(abs.(diff)))
    end

    return max_diff
end

function create_psf_comparison(scalar_psf, spline_psf, xgrid, ygrid, zgrid)
    # Find global intensity range
    local_min, local_max = find_intensity_range(scalar_psf, spline_psf, xgrid, ygrid, zgrid)

    # Find max difference for difference plots
    max_diff = find_max_difference(scalar_psf, spline_psf, xgrid, ygrid, zgrid)

    # Create figure with proper spacing
    fig = Figure(size=(1200, 900), backgroundcolor=:white)  # Increased height for RMSE row

    # Create a proper grid layout with specific sizing
    gl = fig[1, 1] = GridLayout()
    
    # Add title at the top spanning all columns
    title_row = 1
    content_start_row = 2
    gl[title_row, 2:(length(zgrid)+1)] = Label(fig, "PSF Comparison Across Z-Planes", fontsize=24)

    # Add row titles in their own column
    row_label_col = 1
    content_start_col = 2
    gl[content_start_row, row_label_col] = Label(fig, "Original", fontsize=18, rotation=π/2, tellheight=false)
    gl[content_start_row+1, row_label_col] = Label(fig, "Spline", fontsize=18, rotation=π/2, tellheight=false)
    gl[content_start_row+2, row_label_col] = Label(fig, "Difference", fontsize=18, rotation=π/2, tellheight=false)
    gl[content_start_row+3, row_label_col] = Label(fig, "RMSE", fontsize=18, rotation=π/2, tellheight=false)

    # Add column titles
    for (i, z) in enumerate(zgrid)
        gl[title_row+1, content_start_col+i-1] = Label(fig, "z = $(z) μm", fontsize=16)
    end

    # Create plots
    axes = Dict()
    
    # Create two colorbars - one for PSF and one for difference
    colorbar_col = content_start_col + length(zgrid)
    
    # PSF colorbar
    cb_psf = Colorbar(gl[content_start_row:content_start_row+1, colorbar_col], 
                     colormap=:viridis,
                     limits=(local_min, local_max), 
                     labelsize=14, 
                     label="Intensity",
                     height=Relative(0.9))
    
    # Difference colorbar
    cb_diff = Colorbar(gl[content_start_row+2, colorbar_col], 
                      colormap=:balance,
                      limits=(-max_diff, max_diff), 
                      labelsize=14, 
                      label="Difference",
                      height=Relative(0.9))

    for (i, z) in enumerate(zgrid)
        col = content_start_col + i - 1
        
        original = [scalar_psf(x, y, z) for y in ygrid, x in xgrid]
        spline = [spline_psf(x, y, z) for y in ygrid, x in xgrid]
        diff = original - spline

        # Original PSF
        axes[:orig, i] = Axis(gl[content_start_row, col], aspect=DataAspect())
        hidedecorations!(axes[:orig, i])
        hm1 = heatmap!(axes[:orig, i], xgrid, ygrid, original, colormap=:viridis,
            colorrange=(local_min, local_max))

        # Spline PSF
        axes[:spline, i] = Axis(gl[content_start_row+1, col], aspect=DataAspect())
        hidedecorations!(axes[:spline, i])
        hm2 = heatmap!(axes[:spline, i], xgrid, ygrid, spline, colormap=:viridis,
            colorrange=(local_min, local_max))

        # Difference
        axes[:diff, i] = Axis(gl[content_start_row+2, col], aspect=DataAspect())
        hidedecorations!(axes[:diff, i])
        hm3 = heatmap!(axes[:diff, i], xgrid, ygrid, diff, colormap=:balance,
            colorrange=(-max_diff, max_diff))

        # Calculate RMSE
        rmse = sqrt(mean(diff .^ 2))
        
        # Add RMSE values in separate row
        rmse_label = Label(gl[content_start_row+3, col], 
                          "$(round(rmse, digits=5))",
                          fontsize=14,
                          tellwidth=false)

        println("  z = $(z)μm: RMSE = $(rmse)")
    end

    # Set consistent row heights and column widths
    for row in content_start_row:(content_start_row+2)
        rowsize!(gl, row, Relative(0.27))
    end
    rowsize!(gl, content_start_row+3, Relative(0.1))  # Smaller height for RMSE row
    
    for col in content_start_col:(content_start_col+length(zgrid)-1)
        colsize!(gl, col, Relative(1.0/length(zgrid)))
    end
    
    # Apply proper spacing
    rowgap!(gl, 10)
    colgap!(gl, 10)

    return fig
end

# Use the function to create and display the visualization
fig = create_psf_comparison(scalar_psf, spline_psf, xgrid, ygrid, zgrid)
display(fig)

#= PART 2: PERFORMANCE COMPARISON =#

# Create a test grid
test_coords = [(rand() - 0.5, rand() - 0.5, (rand() - 0.5) * 0.6) for _ in 1:1000]

# Benchmark original PSF
println("\nBenchmarking original PSF evaluation...")
orig_time = @benchmark begin
    for (x, y, z) in $test_coords
        $scalar_psf(x, y, z)
    end
end

# Benchmark spline PSF
println("Benchmarking spline PSF evaluation...")
spline_time = @benchmark begin
    for (x, y, z) in $test_coords
        $spline_psf(x, y, z)
    end
end

println("\nPerformance summary:")
println("  Original PSF: $(mean(orig_time).time / 1e6) ms")
println("  Spline PSF: $(mean(spline_time).time / 1e6) ms")
println("  Speedup: $(mean(orig_time).time / mean(spline_time).time)x faster")

