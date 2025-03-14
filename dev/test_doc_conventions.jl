using MicroscopePSFs

println("Testing conventions.md")

# Type Hierarchy
println("\nChecking type hierarchy...")

# Check if base types exist
println("AbstractPSF: ", @isdefined(AbstractPSF))
println("Abstract2DPSF: ", @isdefined(Abstract2DPSF))
println("Abstract3DPSF: ", @isdefined(Abstract3DPSF))

# Check if specific PSF types exist and are subtypes
println("Gaussian2D: ", @isdefined(Gaussian2D))
println("Airy2D: ", @isdefined(Airy2D))
println("Scalar3DPSF: ", @isdefined(Scalar3DPSF))
println("Vector3DPSF: ", @isdefined(Vector3DPSF))
println("SplinePSF: ", @isdefined(SplinePSF))

# Check if 2D PSFs are subtypes of Abstract2DPSF
if @isdefined(Gaussian2D) && @isdefined(Abstract2DPSF)
    println("Gaussian2D <: Abstract2DPSF: ", Gaussian2D <: Abstract2DPSF)
end

if @isdefined(Airy2D) && @isdefined(Abstract2DPSF)
    println("Airy2D <: Abstract2DPSF: ", Airy2D <: Abstract2DPSF)
end

# Check if 3D PSFs are subtypes of Abstract3DPSF
if @isdefined(Scalar3DPSF) && @isdefined(Abstract3DPSF)
    println("Scalar3DPSF <: Abstract3DPSF: ", Scalar3DPSF <: Abstract3DPSF)
end

if @isdefined(Vector3DPSF) && @isdefined(Abstract3DPSF)
    println("Vector3DPSF <: Abstract3DPSF: ", Vector3DPSF <: Abstract3DPSF)
end

println("\nAll checks completed for conventions.md")