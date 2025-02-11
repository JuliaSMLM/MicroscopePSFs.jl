# dev/test_vector3d.jl

using Revise 
using MicroscopePSFs
using CairoMakie
using Printf

# Basic parameters
λ = 0.532        # Wavelength (μm)
na = 1.4         # Numerical aperture
n_medium = 1.33  # Sample medium (water)
n_coverslip = 1.52  # Cover slip (glass)
n_immersion = 1.52  # Immersion medium (oil)

# Create three orthogonal dipoles
dipole_x = DipoleVector(1.0, 0.0, 0.0)
dipole_y = DipoleVector(0.0, 1.0, 0.0)
dipole_z = DipoleVector(0.0, 0.0, 1.0)

# Add some aberration via Zernike coefficients
zc = ZernikeCoefficients(10)  # First 10 terms
zc.phase[4] = 0.5  # Add some defocus
zc.phase[8] = 0.2  # Add some spherical aberration

# Create PSFs
psf_x = Vector3DPSF(na, λ, dipole_x; 
    n_medium=n_medium, 
    n_coverslip=n_coverslip,
    n_immersion=n_immersion,
    base_zernike=zc)

psf_y = Vector3DPSF(na, λ, dipole_y;
    n_medium=n_medium,
    n_coverslip=n_coverslip,
    n_immersion=n_immersion,
    base_zernike=zc)

psf_z = Vector3DPSF(na, λ, dipole_z;
    n_medium=n_medium,
    n_coverslip=n_coverslip,
    n_immersion=n_immersion,
    base_zernike=zc)

# Create random orientation by averaging
psf_random(x, y, z) = (psf_x(x, y, z) + psf_y(x, y, z) + psf_z(x, y, z))/3

# Create visualization grids
x = y = range(-1, 1, 100)
z = range(-2, 2, 21)  # Axial range for XZ view

# Compute intensities
intensity_xy_x = [psf_x(xi, yi, 0.0) for yi in y, xi in x]
intensity_xy_y = [psf_y(xi, yi, 0.0) for yi in y, xi in x]
intensity_xy_z = [psf_z(xi, yi, 0.0) for yi in y, xi in x]
intensity_xy_random = [psf_random(xi, yi, 0.0) for yi in y, xi in x]

# Compute XZ sections
intensity_xz_x = [psf_x(xi, 0.0, zi) for zi in z, xi in x]
intensity_xz_random = [psf_random(xi, 0.0, zi) for zi in z, xi in x]

# Create figure
fig = Figure(size=(1200, 800))

# Helper function for consistent axis formatting
function setup_axis!(ax, title)
    ax.aspect = DataAspect()
    ax.xlabel = "x (μm)"
    ax.title = title
end

# XY slices for different dipole orientations
ax1 = Axis(fig[1, 1]; ylabel="y (μm)")
setup_axis!(ax1, "X dipole")
hm1 = heatmap!(ax1, x, y, intensity_xy_x, colormap=:viridis)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3])
setup_axis!(ax2, "Y dipole")
hm2 = heatmap!(ax2, x, y, intensity_xy_y, colormap=:viridis)
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1]; ylabel="y (μm)")
setup_axis!(ax3, "Z dipole")
hm3 = heatmap!(ax3, x, y, intensity_xy_z, colormap=:viridis)
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3])
setup_axis!(ax4, "Random orientation")
hm4 = heatmap!(ax4, x, y, intensity_xy_random, colormap=:viridis)
Colorbar(fig[2, 4], hm4)

# XZ sections
ax5 = Axis(fig[3, 1]; ylabel="z (μm)")
setup_axis!(ax5, "XZ section - X dipole")
hm5 = heatmap!(ax5, x, z, intensity_xz_x, colormap=:viridis)
Colorbar(fig[3, 2], hm5)

ax6 = Axis(fig[3, 3]; ylabel="z (μm)")
setup_axis!(ax6, "XZ section - Random")
hm6 = heatmap!(ax6, x, z, intensity_xz_random, colormap=:viridis)
Colorbar(fig[3, 4], hm6)

# Add descriptive title
Label(fig[0, :], text="Vector PSF Tests (NA=$na, λ=$λ μm, n_medium=$n_medium)",
    fontsize=20)

# Print validation info
println("\nValidation:")
println("PSF parameters:")
println("  NA = ", na)
println("  λ = ", λ, " μm")
println("  n_medium = ", n_medium)
println("  n_coverslip = ", n_coverslip)
println("  n_immersion = ", n_immersion)

println("\nGrid parameters:")
println("  XY range: [", extrema(x), "] μm")
println("  Z range: [", extrema(z), "] μm")
println("  Grid size: ", length(x), "x", length(y))

println("\nEnergy conservation check:")
total_energy_x = sum(intensity_xy_x) * (x[2]-x[1]) * (y[2]-y[1])
total_energy_random = sum(intensity_xy_random) * (x[2]-x[1]) * (y[2]-y[1])
println("  X dipole total energy: ", @sprintf("%.6f", total_energy_x))
println("  Random dipole total energy: ", @sprintf("%.6f", total_energy_random))

# Save figure
save("dev/vector3d_test.png", fig)

# Return figure handle for display
fig