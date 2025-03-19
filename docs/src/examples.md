# Examples

This page provides practical examples demonstrating key features of MicroscopePSFs.jl.

## 3D Multi-Emitter Simulation

This example shows how to efficiently simulate many 3D emitters using VectorPSF with astigmatism, SplinePSF for acceleration, and support regions for performance optimization:

```jldoctest; output = false
using MicroscopePSFs
using Random

# 1. Create Zernike coefficients for astigmatism
zc = ZernikeCoefficients(15)  # Up to 15th Zernike term
zc.phase[6] = 0.5  # Add astigmatism (OSA indexing + 1, normalized to rms = 1.0)

# 2. Create a rotating dipole VectorPSF with astigmatism
# No explicit dipole orientation means a rotating dipole is used
vector_psf = VectorPSF(
    1.4,            # NA
    0.680,          # Wavelength (μm)
    n_medium=1.33,  # Water
    n_immersion=1.518,
    base_zernike=zc  # Apply the astigmatism
)

# 3. Create faster SplinePSF approximation
# This pre-computes the PSF on a grid for much faster evaluation
x_range = y_range = range(-2.0, 2.0, step=0.1)
z_range = range(-2.0, 2.0, step=0.2)
spline_psf = SplinePSF(vector_psf, x_range, y_range, z_range)

# 4. Create camera (100x100 pixels, 100nm pixel size)
camera = IdealCamera(100, 100, 0.1)  # nx, ny, pixel_size in μm

# 5. Create multiple random 3D emitters
Random.seed!(42)
n_emitters = 100
dipole_z = DipoleVector(0.0, 0.0, 1.0)  # For fixed orientation emitters

# Mix of fixed and rotating dipoles
emitters = []
for i in 1:n_emitters
    x = 10*rand()      # x (μm)
    y = 10*rand()      # y (μm)
    z = 2*(rand()-0.5) # z (μm)
    photons = 1000*rand()  # photons
    
    # Every other emitter gets a z-oriented dipole, others get random orientation
    if i % 2 == 0
        push!(emitters, DipoleEmitter3D(x, y, z, photons, dipole_z))
    else
        # Random dipole orientation
        θ = π*rand()     # polar angle
        φ = 2π*rand()    # azimuthal angle
        dx = sin(θ)*cos(φ)
        dy = sin(θ)*sin(φ)
        dz = cos(θ)
        push!(emitters, DipoleEmitter3D(x, y, z, photons, dx, dy, dz))
    end
end

# 6. Generate camera image with all emitters
# Using support=1.0 dramatically speeds up the calculation
# by only computing each PSF within a 1μm radius
pixels = integrate_pixels(spline_psf, camera, emitters; support=1.0)

# Display result (using CairoMakie in real application)
size(pixels)

# output
(100, 100)
```
