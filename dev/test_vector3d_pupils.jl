# dev/test_vector3d.jl

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
using SMLMData

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index of immersion medium
n_med = 1.33  # Refractive index of medium (water)

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(1:nx, 1:ny, pixel_size)
cam_size = nx * pixel_size

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
# Add astigmatism
zc.phase[6] = 1.0

# Create base pupil function
base_pupil = PupilFunction(na, λ, n, zc)

# Create Vector PSF
psf = Vector3DPupilPSF(na, λ, base_pupil; n_medium=n_med)

# Evaluation coordinates
x = y = range(-1, 1, 100)  # PSF field coordinates
z_planes = [-1.0, -0.5, 0.0, 0.5, 1.0]  # z positions in microns

# Different dipole orientations to test
dipoles = [
    DipoleVector(1.0, 0.0, 0.0),  # x-oriented
    DipoleVector(0.0, 1.0, 0.0),  # y-oriented
    DipoleVector(0.0, 0.0, 1.0)  # z-oriented
]

# Create visualization
fig = Figure(
    size = (2400, 1200),
    backgroundcolor = :white,
    figure_padding = (40, 40, 40, 40)
)

# Helper function for axis setup
function setup_axis!(ax, title)
    ax.title = title
    ax.xlabel = "x index"
    ax.ylabel = "y index"
    ax.aspect = DataAspect()
    ax.yreversed = true  # Put (1,1) at top left
end

# Get all pupils from the vector pupil
vpupil = psf.pupil
pupils = vpupil.pupils
dipole_labels = ["x-dipole", "y-dipole", "z-dipole"]

# Create Ex plots (left half)
for i in 1:3
    # Magnitude plot for Ex
    ax_mag = Axis(fig[i, 1], width=225, height=225)
    setup_axis!(ax_mag, "Ex from $(dipole_labels[i])\nMagnitude")
    hm_mag = heatmap!(ax_mag, 1:length(x), 1:length(y), abs.(pupils[i].field)', colormap=:viridis)
    i == 1 && (Colorbar(fig[i, 2], hm_mag))

    # Phase plot for Ex
    ax_phase = Axis(fig[i, 3], width=225, height=225)
    setup_axis!(ax_phase, "Ex from $(dipole_labels[i])\nPhase")
    hm_phase = heatmap!(ax_phase, 1:length(x), 1:length(y), angle.(pupils[i].field)', 
                        colormap=:twilight, colorrange=(-π, π))
    i == 1 && (Colorbar(fig[i, 4], hm_phase))

    # Magnitude plot for Ey
    ax_mag = Axis(fig[i, 5], width=225, height=225)
    setup_axis!(ax_mag, "Ey from $(dipole_labels[i])\nMagnitude")
    hm_mag = heatmap!(ax_mag, 1:length(x), 1:length(y), abs.(pupils[i+3].field)', colormap=:viridis)
    i == 1 && (Colorbar(fig[i, 6], hm_mag))

    # Phase plot for Ey
    ax_phase = Axis(fig[i, 7], width=225, height=225)
    setup_axis!(ax_phase, "Ey from $(dipole_labels[i])\nPhase")
    hm_phase = heatmap!(ax_phase, 1:length(x), 1:length(y), angle.(pupils[i+3].field)', 
                        colormap=:twilight, colorrange=(-π, π))
    i == 1 && (Colorbar(fig[i, 8], hm_phase))
end

# Add column labels
Label(fig[0, 1:4], "Ex Fields", fontsize=20)
Label(fig[0, 5:8], "Ey Fields", fontsize=20)

# Add row labels
Label(fig[1:1, 0], "x-dipole", rotation=π/2, fontsize=16)
Label(fig[2:2, 0], "y-dipole", rotation=π/2, fontsize=16)
Label(fig[3:3, 0], "z-dipole", rotation=π/2, fontsize=16)

# Add title
Label(fig[0, :], "Vector Pupil Fields\n(NA=$na, λ=$λ μm, n=$n, n_medium=$n_med)", 
      fontsize=20, padding=(0,0,20,0))

# Adjust layout
rowgap!(fig.layout, 20)
colgap!(fig.layout, 20)

# Save figure
save("dev/vector_pupils.png", fig)

# Print parameters
println("\nPSF parameters:")
println("NA = $na, λ = $λ μm")
println("n_immersion = $n, n_medium = $n_med")
println("Astigmatism: 1 wave")

fig