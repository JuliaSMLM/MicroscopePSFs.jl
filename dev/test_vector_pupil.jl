using Pkg 
Pkg.activate("dev")
using Revise
using MicroscopePSFs
using CairoMakie
using Printf

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n_medium = 1.33   # Sample medium refractive index (water)
n_coverslip = 1.52  # Coverslip refractive index (glass)
n_immersion = 1.52  # Immersion medium refractive index (oil)
z_stage = 0.0  # Distance the sample stage was moved away from the nominal focal plane at the coverslip

# Create dipole orientations
dipole_x = DipoleVector(1.0, 0.0, 0.0)  # x-oriented
dipole_y = DipoleVector(0.0, 1.0, 0.0)  # y-oriented
dipole_z = DipoleVector(0.0, 0.0, 1.0)  # z-oriented
dipoles = [dipole_x, dipole_y, dipole_z]
dipole_labels = ["x-dipole", "y-dipole", "z-dipole"]

# Grid size for pupil visualization
grid_size = 128

# Create Vector3DPSF for each dipole orientation
psfs = [
    Vector3DPSF(na, λ, dipole; 
                n_medium=n_medium, 
                n_coverslip=n_coverslip, 
                n_immersion=n_immersion,
                z_stage=z_stage,
                grid_size=grid_size)
    for dipole in dipoles
]

# Create figure with 3 rows (dipoles) and 4 columns (Ex mag, Ex phase, Ey mag, Ey phase)
fig = Figure(size=(1600, 1200), 
             backgroundcolor = :white,
             figure_padding = 30)

# Normalize colormap ranges for consistent visualization
ex_max = maximum([maximum(abs.(psf.vector_pupils.Ex.field)) for psf in psfs])
ey_max = maximum([maximum(abs.(psf.vector_pupils.Ey.field)) for psf in psfs])

# Column titles
column_titles = ["Ex Magnitude", "Ex Phase", "Ey Magnitude", "Ey Phase"]
for (i, title) in enumerate(column_titles)
    Label(fig[0, i], title, fontsize=16)
end

# Function to hide axis decorations
function clean_axis!(ax)
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xlabelvisible = false
    ax.ylabelvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    # Hide spines
    ax.leftspinevisible = false
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.bottomspinevisible = false
    return ax
end

# Create heatmaps for each dipole orientation and field component
for (i, psf) in enumerate(psfs)
    # Extract pupil fields
    ex_field = psf.vector_pupils.Ex.field
    ey_field = psf.vector_pupils.Ey.field
    
    # Create coordinate grid
    ρ_max = 1.0  # Normalized pupil radius
    xs = ys = range(-ρ_max, ρ_max, grid_size)
    
    # Create a circular mask for the pupil
    mask = [x^2 + y^2 <= ρ_max^2 for y in ys, x in xs]
    
    # Apply mask to fields - set values outside pupil to NaN for nicer visualization
    ex_mag = abs.(ex_field)
    ex_phase = angle.(ex_field)
    ey_mag = abs.(ey_field)
    ey_phase = angle.(ey_field)
    
    for j in 1:grid_size, k in 1:grid_size
        if !mask[j, k]
            ex_mag[j, k] = NaN
            ex_phase[j, k] = NaN
            ey_mag[j, k] = NaN
            ey_phase[j, k] = NaN
        end
    end
    
    # Create heatmaps with manually hidden decorations
    ax_ex_mag = Axis(fig[i, 1], aspect=DataAspect())
    clean_axis!(ax_ex_mag)
    heatmap!(ax_ex_mag, xs, ys, ex_mag, colormap=:viridis, colorrange=(0, ex_max))
    
    ax_ex_phase = Axis(fig[i, 2], aspect=DataAspect())
    clean_axis!(ax_ex_phase)
    heatmap!(ax_ex_phase, xs, ys, ex_phase, colormap=:twilight, colorrange=(-π, π))
    
    ax_ey_mag = Axis(fig[i, 3], aspect=DataAspect())
    clean_axis!(ax_ey_mag)
    heatmap!(ax_ey_mag, xs, ys, ey_mag, colormap=:viridis, colorrange=(0, ey_max))
    
    ax_ey_phase = Axis(fig[i, 4], aspect=DataAspect())
    clean_axis!(ax_ey_phase)
    heatmap!(ax_ey_phase, xs, ys, ey_phase, colormap=:twilight, colorrange=(-π, π))
    
    # Add row label for dipole orientation
    Label(fig[i, 0], dipole_labels[i], rotation=π/2, padding=(0, 10), fontsize=16)
end

# Add overall title
title_text = @sprintf("Pupil Functions for Different Dipole Orientations\n(NA=%.1f, λ=%.3f μm, n_medium=%.2f)", 
                     na, λ, n_medium)
Label(fig[0, 1:4], title_text, fontsize=20, font=:bold, padding=(0, 0, 20, 0))

# Adjust layout
rowgap!(fig.layout, 20)
colgap!(fig.layout, 20)

# Save figure
save("dipole_pupils_clean.png", fig, px_per_unit=2)  # Higher resolution

# Display information
println("Pupil Parameters:")
println("NA = $na, λ = $λ μm")
println("Medium index = $n_medium")
println("Coverslip index = $n_coverslip")
println("Immersion index = $n_immersion")
println("Z stage position = $z_stage μm")
println("Grid size = $grid_size × $grid_size")

# Return figure for display
fig