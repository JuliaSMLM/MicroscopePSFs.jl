# check_type_stability.jl
# Script to check type stability for all PSF types
using Pkg 
Pkg.activate("dev")

using MicroscopePSFs
using InteractiveUtils
using Enzyme

# Function to run type checks
function run_check(psf_name, psf, emitter_name, emitter)
    println("="^80)
    println("Testing: $psf_name with $emitter_name")
    println("-"^80)
    
    # Create camera
    nx = ny = 20  # 20x20 camera
    pixel_size = 0.1  # 100 nm pixels
    camera = IdealCamera(nx, ny, pixel_size)
    
    # Call the function and print warnings
    println("## @code_warntype for integrate_pixels:")
    @code_warntype integrate_pixels(psf, camera, emitter)
    
    println("\n## @code_warntype for integrate_pixels with support:")
    @code_warntype integrate_pixels(psf, camera, emitter, support=1.0)
    
    println("\n## @code_warntype for Enzyme autodiff:")
    @code_warntype autodiff(Forward, x -> sum(integrate_pixels(psf, camera, 
        typeof(emitter)(x, emitter.y, emitter.z, emitter.photons))), emitter.x)
    
    println("="^80)
    println()
end

# Create emitters
emitter2D = Emitter2D(1.0, 1.0, 1000.0)
emitter3D = Emitter3D(1.0, 1.0, 0.5, 1000.0)
dipole_emitter = DipoleEmitter3D(1.0, 1.0, 0.5, 1000.0, 1.0, 0.0, 0.0)

#=====================================
  2D PSF tests
=====================================#

# GaussianPSF
gaussian_psf = GaussianPSF(0.15)  # sigma = 150nm
run_check("GaussianPSF", gaussian_psf, "Emitter2D", emitter2D)
run_check("GaussianPSF", gaussian_psf, "Emitter3D", emitter3D)

# AiryPSF
airy_psf = AiryPSF(1.4, 0.532)  # NA=1.4, λ=532nm
run_check("AiryPSF", airy_psf, "Emitter2D", emitter2D)
run_check("AiryPSF", airy_psf, "Emitter3D", emitter3D)

#=====================================
  3D PSF tests
=====================================#

# ScalarPSF
scalar_psf = ScalarPSF(1.4, 0.532, 1.52)  # NA=1.4, λ=532nm, n=1.52
run_check("ScalarPSF", scalar_psf, "Emitter3D", emitter3D)

# VectorPSF
vector_psf = VectorPSF(1.4, 0.532, DipoleVector(1.0, 0.0, 0.0))
run_check("VectorPSF", vector_psf, "DipoleEmitter3D", dipole_emitter)

# SplinePSF from a ScalarPSF
x_range = y_range = range(-1.0, 1.0, length=21)
z_range = range(-0.5, 0.5, length=11)
spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)
run_check("SplinePSF", spline_psf, "Emitter3D", emitter3D)

#=====================================
  Batch emitter tests
=====================================#

# Test with multiple emitters
emitters = [
    Emitter3D(0.8, 0.8, 0.1, 1000.0),
    Emitter3D(1.2, 1.2, 0.2, 1000.0)
]

println("="^80)
println("Testing: ScalarPSF with multiple emitters")
println("-"^80)
@code_warntype integrate_pixels(scalar_psf, camera, emitters)
println("\nWith support:")
@code_warntype integrate_pixels(scalar_psf, camera, emitters, support=1.0)
println("="^80)