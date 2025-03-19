# dev/test_scalar3d.jl
using Pkg 
Pkg.activate("dev")

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
# using Zygote 
using Enzyme


# Microscope parameters
λ = 0.532  # Green light wavelength in microns
na = 1.4   # Numerical aperture
n = 1.52   # Refractive index

# Camera setup
pixel_size = 0.1  # 100 nm pixels
nx = 20           # Width
ny = 10           # Height
camera = IdealCamera(1:nx, 1:ny, pixel_size)

# Create PSF with astigmatism
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.0  # Add phase shift for visualization
psf = ScalarPSF(na, λ, n, coeffs=zc)

# Calculate the derivative with respect to x
# dpsf_dx = x -> Zygote.gradient(x -> psf(x, y_value, z_value), x)[1]

# Calc derivative image 
x_range = y_range = range(-0.5, 0.5, 20)  # PSF field coordinates
# Only compute the derivative with respect to x

# Only compute the derivative with respect to x


dpsf = Enzyme.autodiff(Enzyme.Reverse,  x -> psf(x, 0.0, 0.0), Active, Active(1.0))
dx_values = [Enzyme.autodiff(Enzyme.Reverse, x -> psf(x, y, 0.0), Active, Active(x))[1] for x in x_range, y in y_range]

# integrate pixels
dx_image = [Enzyme.autodiff(Enzyme.Reverse, x -> integrate_pixels(psf, camera, Emitter3D(x, y, 0.0,1000.0)), Active, Active(x))[1] for x in x_range, y in y_range]


dx_values = [gradient(x -> psf(x, y, 0.0), x)[1] for x in x_range, y in y_range]
dy_values = [gradient(y -> psf(x, y, 0.0), y)[1] for x in x_range, y in y_range]
dz_values = [gradient(z -> psf(x, y, z), 1.0)[1] for x in x_range, y in y_range]
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

# Now try with integrating the PSF
dx_image = [gradient(x -> integrate_pixels(psf, camera, Emitter3D(x, y, 0.0,1000.0)), x)[1] for x in x_range, y in y_range]
dy_image = [gradient(y -> integrate_pixels(psf, camera, Emitter3D(x, y, 0.0,1000.0)), y)[1] for x in x_range, y in y_range]
dz_image = [gradient(z -> integrate_pixels(psf, camera, Emitter3D(x, y, z,1000.0)), 1.0)[1] for x in x_range, y in y_range]

image = [integrate_psf(psf, x, y, 0.0) for x in x_range, y in y_range]

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









