using Pkg 
Pkg.activate("dev")
using Revise 
using MicroscopePSFs
using CairoMakie
using LinearAlgebra
using Statistics
using Zygote
using BenchmarkTools

#= PART 1: BASIC USAGE OF SPLINEPSF =#


# Create a Scalar3D PSF
scalar_psf = Scalar3DPSF(1.4, 0.532, 1.518)

# Sample to create a SplinePSF
x_range = y_range = range(-2.0, 2.0, length=41)  # 41×41 lateral grid, 4μm range
z_range = range(-1.0, 1.0, length=21)            # 21 z-planes, 2μm range

println("Type of x_range: ", typeof(x_range))
println("Type of y_range: ", typeof(y_range))
println("Type of z_range: ", typeof(z_range))

# Update the SplinePSF creation to use the positional constructor
println("Creating SplinePSF from Scalar3DPSF...")
spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)

println("Original PSF type: $(typeof(scalar_psf))")
println("Spline PSF type: $(typeof(spline_psf))")

# Create grids for visualization
xgrid = ygrid = range(-1.0, 1.0, length=101)  # Finer grid for visualization
zgrid = [0.0, 0.2, 0.4, 0.6]  # Few z-planes to visualize

# Compare PSFs at different z-planes
println("Comparing PSFs at different z-planes...")
figs = []  # To store the figures
for z in zgrid
    # Calculate PSF intensities
    original = [scalar_psf(x, y, z) for y in ygrid, x in xgrid]
    spline   = [spline_psf(x, y, z) for y in ygrid, x in xgrid]
    
    # Create a new figure for this z-plane.
    # Here we allocate two rows: the first for a global title and the second for the two side-by-side axes.
    fig = Figure(resolution = (800, 400))
    fig[1, :] = Label(fig, "z = $(z) μm", fontsize=20, tellwidth=false)
    
    ax1 = Axis(fig[2, 1], title = "Original", aspect = DataAspect())
    ax2 = Axis(fig[2, 2], title = "Spline", aspect = DataAspect())
    
    heatmap!(ax1, xgrid, ygrid, original)
    heatmap!(ax2, xgrid, ygrid, spline)
    
    push!(figs, fig)
    
    # Calculate difference and RMSE
    diff = original - spline
    rmse = sqrt(mean(diff.^2))
    println("  z = $(z)μm: RMSE = $(rmse)")
end

# Optionally, display one of the figures (when running interactively)
# For example, to display the first figure:
# display(figs[1])

#= PART 2: PERFORMANCE COMPARISON =#


# Create a test grid
test_coords = [(rand()-0.5, rand()-0.5, (rand()-0.5)*0.6) for _ in 1:1000]

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

#= PART 3: AUTOMATIC DIFFERENTIATION =#


# Simple loss function to minimize
function psf_loss(params, psf)
    x, y, z = params
    return -psf(x, y, z)  # Negative because we want to maximize intensity
end

# Initial guess
initial_params = [0.1, 0.2, 0.3]

# Calculate gradient with Zygote
println("\nCalculating gradient with Zygote...")
loss_orig(p) = psf_loss(p, scalar_psf)
loss_spline(p) = psf_loss(p, spline_psf)

grad_orig = Zygote.gradient(loss_orig, initial_params)[1]
grad_spline = Zygote.gradient(loss_spline, initial_params)[1]

println("Gradients at $(initial_params):")
println("  Original PSF: $(grad_orig)")
println("  Spline PSF: $(grad_spline)")
println("  Cosine similarity: $(dot(grad_orig, grad_spline) / (norm(grad_orig) * norm(grad_spline)))")

# Gradient descent
function gradient_descent(loss_fn, initial_params; steps=10, learning_rate=0.1)
    params = copy(initial_params)
    trajectory = [copy(params)]
    
    for i in 1:steps
        grad = Zygote.gradient(loss_fn, params)[1]
        params .-= learning_rate .* grad
        push!(trajectory, copy(params))
    end
    
    return trajectory
end

println("\nRunning gradient descent optimization...")
traj_orig = gradient_descent(loss_orig, initial_params)
traj_spline = gradient_descent(loss_spline, initial_params)

println("Final positions:")
println("  Original PSF: $(traj_orig[end])")
println("  Spline PSF: $(traj_spline[end])")
println("  Distance: $(norm(traj_orig[end] - traj_spline[end]))")

println("\nFinal intensities:")
println("  Original PSF: $(-loss_orig(traj_orig[end]))")
println("  Spline PSF: $(-loss_spline(traj_spline[end]))")
