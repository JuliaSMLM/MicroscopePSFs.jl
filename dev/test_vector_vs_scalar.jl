# defocus_comparison.jl
using Pkg
Pkg.activate(".")

using MicroscopePSFs
using CairoMakie
using Printf
using SMLMData

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n_medium = 1.52   # Sample medium refractive index (water)
n_coverslip = 1.52  # Coverslip refractive index (glass)
n_immersion = 1.52  # Immersion medium refractive index (oil)
focal_z = 1.0  # Nominal focal plane position

# Create dipole orientation (x-oriented)
dipole = DipoleVector(1.0, 1.0, 0.0)

# Create Zernike coefficients with astigmatism
zc = ZernikeCoefficients(15)
zc.phase[6] = 1.0  # Add vertical astigmatism

# Create Vector3DPSF with x-oriented dipole
vector_psf = Vector3DPSF(
    na, λ, dipole; 
    base_zernike=zc,
    n_medium=n_medium, 
    n_coverslip=n_coverslip, 
    n_immersion=n_immersion,
    focal_z=focal_z
)

# Create Scalar3DPSF with the same aberration
scalar_psf = Scalar3DPSF(na, λ, n_medium; coeffs=zc)

# Evaluation coordinates
x = y = range(-1, 1, 100)  # PSF field coordinates
z_planes = collect(-1.0:0.5:1.0)  # z positions in microns

# Create figure with a nice layout
fig = Figure(
    size = (1500, 800),
    backgroundcolor = RGBf(0.98, 0.98, 0.98),
    figure_padding = (40, 40, 40, 40)
)

# Create main GridLayouts
g_scalar = fig[1, 1] = GridLayout()
g_vector = fig[2, 1] = GridLayout()
g_slice = fig[1:2, 2] = GridLayout()

# Helper function for axis setup
function setup_axis!(ax, title)
    ax.title = title
    ax.xlabel = "x (μm)"
    ax.ylabel = "y (μm)"
    ax.aspect = DataAspect()
end

# Create axes for scalar PSFs
scalar_axes = [Axis(g_scalar[1, i]) for i in 1:length(z_planes)]
for (i, (ax, z)) in enumerate(zip(scalar_axes, z_planes))
    setup_axis!(ax, @sprintf("z = %.1f μm", z))
    
    # Calculate scalar PSF intensity
    println("Calculating scalar PSF at z = $z")
    scalar_intensity = [scalar_psf(xi, yi, z) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, scalar_intensity, colormap=:viridis)
    
    # Add colorbar for last plot
    if i == length(z_planes)
        Colorbar(g_scalar[1, i+1], hm, width=15, label="Intensity")
    end
end

# Create axes for vector PSFs
vector_axes = [Axis(g_vector[1, i]) for i in 1:length(z_planes)]
for (i, (ax, z)) in enumerate(zip(vector_axes, z_planes))
    setup_axis!(ax, @sprintf("z = %.1f μm", z))
    
    # Calculate vector PSF intensity
    vector_intensity = [vector_psf(xi, yi, z) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, vector_intensity, colormap=:viridis)
    
    # Add colorbar for last plot
    if i == length(z_planes)
        Colorbar(g_vector[1, i+1], hm, width=15, label="Intensity")
    end
end

# Create cross-section plots
ax_xslice = Axis(g_slice[1, 1], 
    title = "PSF Cross-Section at y=0",
    xlabel = "x (μm)",
    ylabel = "Intensity")

ax_zslice = Axis(g_slice[2, 1], 
    title = "PSF Axial Profile at (x,y)=(0,0)",
    xlabel = "z (μm)",
    ylabel = "Intensity")

# Calculate cross-section data
x_profile = [scalar_psf(xi, 0.0, 0.0) for xi in x]
v_x_profile = [vector_psf(xi, 0.0, 0.0) for xi in x]

z_fine = range(-1.0, 1.0, 100)
z_profile = [scalar_psf(0.0, 0.0, zi) for zi in z_fine]
v_z_profile = [vector_psf(0.0, 0.0, zi) for zi in z_fine]

# Plot cross-sections
lines!(ax_xslice, x, x_profile, label="Scalar PSF", linewidth=3)
lines!(ax_xslice, x, v_x_profile, label="Vector PSF", linewidth=3, linestyle=:dash)
axislegend(ax_xslice, position=:rt)

lines!(ax_zslice, z_fine, z_profile, label="Scalar PSF", linewidth=3)
lines!(ax_zslice, z_fine, v_z_profile, label="Vector PSF", linewidth=3, linestyle=:dash)
axislegend(ax_zslice, position=:rt)

# Add panel labels
Label(g_scalar[1, 1, TopLeft()], "A",
    fontsize = 24, font = :bold,
    padding = (0, 5, 5, 0), halign = :right)

Label(g_vector[1, 1, TopLeft()], "B",
    fontsize = 24, font = :bold,
    padding = (0, 5, 5, 0), halign = :right)

Label(g_slice[1, 1, TopLeft()], "C",
    fontsize = 24, font = :bold,
    padding = (0, 5, 5, 0), halign = :right)

# Add row labels
Label(g_scalar[1, 0], "Scalar PSF", 
    rotation = π/2, fontsize = 16, 
    padding = (0, 20, 0, 0))

Label(g_vector[1, 0], "Vector PSF", 
    rotation = π/2, fontsize = 16, 
    padding = (0, 20, 0, 0))

# Add main title
title_text = @sprintf("Comparison of Scalar and Vector PSF Models with Astigmatism\n(NA=%.1f, λ=%.3f μm, n_medium=%.2f, dipole=x-oriented)", 
                      na, λ, n_medium)
Label(fig[0, :], title_text, 
    fontsize = 20, font = :bold, 
    padding = (0, 0, 20, 0))

# Adjust layout
colgap!(g_scalar, 10)
colgap!(g_vector, 10)
rowgap!(fig.layout, 30)
colgap!(fig.layout, 30)

# Make the cross-section column less wide
colsize!(fig.layout, 2, Auto(0.4))

# Save figure
save("dev/test_output/scalar_vector_psf_comparison.png", fig)

# Display details
println("PSF Parameters:")
println("NA = $na, λ = $λ μm")
println("Medium index = $n_medium")
println("Coverslip index = $n_coverslip")
println("Immersion index = $n_immersion")
println("Focal plane = $focal_z μm")
println("Astigmatism: 1 wave at 0 degrees")
println("Dipole orientation: ($(dipole.px), $(dipole.py), $(dipole.pz))")
println("Defocus range: $(z_planes[1]) to $(z_planes[end]) μm")

fig