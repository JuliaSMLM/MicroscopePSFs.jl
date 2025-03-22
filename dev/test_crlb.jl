using Pkg 
Pkg.activate("dev")
using Revise
using MicroscopePSFs
using CairoMakie
using Enzyme 

## CRLB Calculation 
using LinearAlgebra 

function crlb_poisson(psf, camera, emitter)
    # Calculate the PSF and its derivatives using Enzyme with the lambda wrapped in Const
    dx_values = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
        Const(x -> integrate_pixels(psf, camera, Emitter3D(x, emitter.y, emitter.z, emitter.photons))),
        emitter.x)[1]
    dy_values = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
        Const(y -> integrate_pixels(psf, camera, Emitter3D(emitter.x, y, emitter.z, emitter.photons))),
        emitter.y)[1]
    dz_values = Enzyme.jacobian(set_runtime_activity(Enzyme.Forward),
        Const(z -> integrate_pixels(psf, camera, Emitter3D(emitter.x, emitter.y, z, emitter.photons))),
        emitter.z)[1]
    
    image = integrate_pixels(psf, camera, Emitter3D(emitter.x, emitter.y, emitter.z, emitter.photons)) .+ eps()
    
    # Calculate the FI matrix
    fi_matrix = zeros(3, 3)
    fi_matrix[1, 1] = sum(dx_values .^ 2 ./ image)
    fi_matrix[2, 1] = fi_matrix[1, 2] = sum(dx_values .* dy_values ./ image)
    fi_matrix[3, 1] = fi_matrix[1, 3] = sum(dx_values .* dz_values ./ image)
    fi_matrix[2, 2] = sum(dy_values .^ 2 ./ image)
    fi_matrix[3, 2] = fi_matrix[2, 3] = sum(dy_values .* dz_values ./ image)
    fi_matrix[3, 3] = sum(dz_values .^ 2 ./ image)

    # Calculate the CRLB
    crlb_matrix = inv(fi_matrix)  # Inverse of the FI matrix
    
    return diag(crlb_matrix)  # Return the diagonal elements (CRLB for x, y, z)
end

# Setup camera
nx, ny, pixel_size = 20, 20, 0.1  # 100 nm pixels
camera = IdealCamera(nx, ny, pixel_size)

# Setup psf 
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5  # Add vertical astigmatism
psf = ScalarPSF(1.2, 0.6, 1.33; zernike_coeffs=zc)

# Setup emitter
x_emitter = y_emitter = nx / 2 * pixel_size
z_emitter, photons = 0.0, 1000.0

# Loop over z positions 
z_range = -1.0:0.05:1.0  # Z stage positions in microns
σ_x = Vector{Float64}(undef, length(z_range))
σ_y = Vector{Float64}(undef, length(z_range))
σ_z = Vector{Float64}(undef, length(z_range))

for (i, z) in enumerate(z_range)
    # Create a new emitter with the current z position
    emitter = Emitter3D(x_emitter, y_emitter, z, photons)
    
    # Calculate the CRLB for the current emitter position
    crlb_values = crlb_poisson(psf, camera, emitter)
    
    σ_x[i] = sqrt(crlb_values[1])  
    σ_y[i] = sqrt(crlb_values[2])
    σ_z[i] = sqrt(crlb_values[3])
end

# Plot with CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], title = "CRLB", xlabel="z (μm)", ylabel="σ (μm)")
lines!(ax, z_range, σ_x, label="σ_x", color=:red)
lines!(ax, z_range, σ_y, label="σ_y", color=:blue)
lines!(ax, z_range, σ_z, label="σ_z", color=:green)
axislegend(ax)
display(fig)
save("crlb.png", fig)
