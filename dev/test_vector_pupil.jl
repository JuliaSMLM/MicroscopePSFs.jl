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

# Create VectorPSF for each dipole orientation
psfs = [
    VectorPSF(na, λ, dipole; 
                n_medium=n_medium, 
                n_coverslip=n_coverslip, 
                n_immersion=n_immersion,
                z_stage=z_stage,
                grid_size=grid_size)
    for dipole in dipoles
]

# Create figure with 3 rows (dipoles) and 4 columns (Ex mag, Ex phase, Ey mag, Ey phase)
fig = Figure(size=(800, 600))

# Define column headers
column_titles = ["Ex Magnitude", "Ex Phase", "Ey Magnitude", "Ey Phase"]

# Add row headers as separate labels
for i in 1:3
    Label(fig[i, 1, Left()], dipole_labels[i], rotation = π/2, padding = (0, 15, 0, 0), fontsize=20)
end

# Create heatmaps for each dipole orientation and field component
for (i, psf) in enumerate(psfs)
    # Extract pupil fields - importantly, direct array access to preserve orientation
    ex_field = psf.vector_pupils[1].Ex.field
    ey_field = psf.vector_pupils[1].Ey.field
    
    # Apply mask to fields - set values outside pupil to NaN for nicer visualization
    ex_mag = abs.(ex_field)
    ex_phase = angle.(ex_field)
    ey_mag = abs.(ey_field)
    ey_phase = angle.(ey_field)
    
    # Create heatmaps with independent scaling for magnitude
    for (j, (data, cmap, crange)) in enumerate([
        (ex_mag, :viridis, (0, maximum(filter(!isnan, ex_mag)))),
        (ex_phase, :twilight, (-π, π)),
        (ey_mag, :viridis, (0, maximum(filter(!isnan, ey_mag)))),
        (ey_phase, :twilight, (-π, π))
    ])
        ax = Axis(fig[i, j], aspect=DataAspect())
        
        # Add column headers to first row only
        if i == 1
            ax.title = column_titles[j]
        end
        
        # Create the heatmap
        heatmap!(ax, data, colormap=cmap, colorrange=crange)
        
        # Hide all decorations
        hidedecorations!(ax)
        hidespines!(ax)
    end
end

# Add overall title
title_text = @sprintf("Pupil Functions for Different Dipole Orientations (NA=%.1f, λ=%.3f μm, n_medium=%.2f)", 
                     na, λ, n_medium)
Label(fig[0, 1:4], title_text, fontsize=16, font=:bold)

# Save figure
# Ensure test_output directory exists
if !isdir("dev/test_output")
    mkdir("dev/test_output")
end
save("dev/test_output/dipole_pupils.png", fig, px_per_unit=2)

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