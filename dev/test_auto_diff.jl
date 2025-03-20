# dev/test_scalar3d.jl
using Pkg 
Pkg.activate("dev")

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
using Enzyme

# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(nx, ny, pixel_size)

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5  # Add vertical astigmatism
psf = ScalarPSF(na, λ, n, coeffs=zc)
psf = VectorPSF(na, λ; base_zernike=zc)  

# Calculate the derivative with respect to x
# dpsf_dx = x -> Zygote.gradient(x -> psf(x, y_value, z_value), x)[1]

# Calc derivative image 
x_range = y_range = range(-0.5, 0.5, 20)  # PSF field coordinates

# Calculate the PSF and its derivatives using Enzyme
dx_values = [Enzyme.autodiff(Enzyme.Reverse, x -> psf(x, y, 0.0), Active, Active(x))[1][1] for x in x_range, y in y_range]
dy_values = [Enzyme.autodiff(Enzyme.Reverse, y -> psf(x, y, 0.0), Active, Active(y))[1][1] for x in x_range, y in y_range]
dz_values = [Enzyme.autodiff(Enzyme.Reverse, z -> psf(x, y, z), Active, Active(0.0))[1][1] for x in x_range, y in y_range]
psf_arr = [psf(x, y, 0.0) for x in x_range, y in y_range]

# Plot with CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], title = "PSF", yreversed=true, aspect=DataAspect())
heatmap!(ax, psf_arr')
ax2 = Axis(fig[1, 2], title = "dPSF/dx", yreversed=true, aspect=DataAspect())
heatmap!(ax2, dx_values')
ax3 = Axis(fig[1, 3], title = "dPSF/dy", yreversed=true, aspect=DataAspect())
heatmap!(ax3, dy_values')
ax4 = Axis(fig[1, 4], title = "dPSF/dz", yreversed=true, aspect=DataAspect())
heatmap!(ax4, dz_values')
display(fig)

# Now try with integrating the PSF. Need to turn off threading for Zygote to work
x_emitter = nx / 2 * pixel_size
y_emitter = ny / 2 * pixel_size
z_emitter = 0.5
photons = 1000.0 

dx_image = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
    x -> integrate_pixels(psf, camera, Emitter3D(x, y_emitter, z_emitter, photons)), 
    x_emitter)[1]
dy_image = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
    y -> integrate_pixels(psf, camera, Emitter3D(x_emitter, y, z_emitter, photons)), 
    y_emitter)[1]
dz_image = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
    z -> integrate_pixels(psf, camera, Emitter3D(x_emitter, y_emitter, z, photons)), 
    z_emitter)[1]

image = integrate_pixels(psf, camera, Emitter3D(x_emitter, y_emitter, z_emitter, photons))

# Plot with CairoMakie
fig2 = Figure()
ax5 = Axis(fig2[1, 1], title = "Integrated PSF", yreversed=true, aspect=DataAspect())
heatmap!(ax5, image')
ax6 = Axis(fig2[1, 2], title = "dPSF/dx", yreversed=true, aspect=DataAspect())
heatmap!(ax6, dx_image')
ax7 = Axis(fig2[1, 3], title = "dPSF/dy", yreversed=true, aspect=DataAspect())
heatmap!(ax7, dy_image')
ax8 = Axis(fig2[1, 4], title = "dPSF/dz", yreversed=true, aspect=DataAspect())
heatmap!(ax8, dz_image')
display(fig2)









