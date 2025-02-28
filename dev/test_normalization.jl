using Pkg
Pkg.activate("dev")
using Revise
using MicroscopePSFs
using SMLMData


# Set up a camera with large enough FOV to capture most of the PSF
# Using 100nm pixels over a 5x5 μm area (50x50 pixels)
pixel_size = 0.1  # μm
nx, ny = 100, 100
x_edges = range(0, nx * pixel_size, nx + 1)
y_edges = range(0, ny * pixel_size, ny + 1)
camera = IdealCamera(collect(x_edges), collect(y_edges))

# Create emitters at the center with 10 photons
center_x, center_y = nx * pixel_size / 2, ny * pixel_size / 2
emitter2d = Emitter2D(center_x, center_y, 10.0)
emitter3d_z0 = Emitter3D(center_x, center_y, 0.0, 10.0)
emitter3d_z1 = Emitter3D(center_x, center_y, 1.0, 10.0)  # 1 μm defocus

# Create dipole for Vector3DPSF
dipole = DipoleVector(1.0, 0.0, 0.0)  # aligned with x axis
emitter_dipole_z0 = DipoleEmitter3D(center_x, center_y, 0.0, 10.0, 1.0, 0.0, 0.0)
emitter_dipole_z1 = DipoleEmitter3D(center_x, center_y, 2.0, 10.0, 1.0, 0.0, 0.0)

println("Testing PSF normalization with 10-photon emitters")
println("=================================================")

# Test 2D PSFs
gaussian2d = Gaussian2D(0.15)  # σ = 150nm
airy2d = Airy2D(1.4, 0.532)    # NA = 1.4, λ = 532nm

pixels_gaussian = integrate_pixels(gaussian2d, camera, emitter2d)
pixels_airy = integrate_pixels(airy2d, camera, emitter2d)

println("2D PSFs:")
println("  Gaussian2D: Sum = $(sum(pixels_gaussian))")
println("  Airy2D: Sum = $(sum(pixels_airy))")

# Test 3D PSFs
scalar3d = Scalar3DPSF(1.4, 0.532, 1.518)  # NA = 1.4, λ = 532nm, n = 1.518
vector3d = Vector3DPSF(1.4, 0.532, dipole)
spline3d = SplinePSF(scalar3d)  # Create spline from scalar3d

# At z=0
pixels_scalar_z0 = integrate_pixels(scalar3d, camera, emitter3d_z0)
pixels_vector_z0 = integrate_pixels(vector3d, camera, emitter_dipole_z0)
pixels_spline_z0 = integrate_pixels(spline3d, camera, emitter3d_z0)

println("\n3D PSFs at z=0:")
println("  Scalar3DPSF: Sum = $(sum(pixels_scalar_z0))")
println("  Vector3DPSF: Sum = $(sum(pixels_vector_z0))")
println("  SplinePSF: Sum = $(sum(pixels_spline_z0))")

# At z=1 μm
pixels_scalar_z1 = integrate_pixels(scalar3d, camera, emitter3d_z1)
pixels_vector_z1 = integrate_pixels(vector3d, camera, emitter_dipole_z1)
pixels_spline_z1 = integrate_pixels(spline3d, camera, emitter3d_z1)

println("\n3D PSFs at z=1 μm:")
println("  Scalar3DPSF: Sum = $(sum(pixels_scalar_z1))")
println("  Vector3DPSF: Sum = $(sum(pixels_vector_z1))")
println("  SplinePSF: Sum = $(sum(pixels_spline_z1))")

function test_depth_dependence(psf::Vector3DPSF, max_depth::Real=5.0)
    depths = range(0, max_depth, length=20)
    intensities = [sum(psf(0.0, 0.0, z)) for z in depths]
    
    # Should show decreasing intensity with depth
    for (z, I) in zip(depths, intensities)
        println("z = $(round(z, digits=2)) μm: I = $(round(I/intensities[1]*100, digits=1))%")
    end
end
test_depth_dependence(vector3d, 5.0)