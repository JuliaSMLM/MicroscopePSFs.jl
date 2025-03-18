#!/usr/bin/env julia
# dev/test_all.jl - Run all test files in the dev directory

using Pkg
Pkg.activate("dev")

using Revise
using MicroscopePSFs
using CairoMakie
using Printf
using BenchmarkTools

# Ensure test_output directory exists
if !isdir("dev/test_output")
    mkdir("dev/test_output")
end

println("== Running all MicroscopePSFs dev tests ==")

# Helper function to run a test file and catch any errors
function run_test_file(filename)
    println("\n\n")
    println("="^50)
    println("Running $filename")
    println("="^50)
    
    try
        include(filename)
        println("✓ $filename completed successfully")
        return true
    catch e
        println("✗ Error in $filename:")
        println(e)
        return false
    end
end

# List of all test files to run
test_files = [
    "test_gaussian2d.jl",      # Start with 2D PSFs
    "test_airy2d.jl",
    "test_Scalar3D.jl",        # Then 3D PSFs
    "test_vector3d.jl",
    "test_vector_pupil.jl",    # Supporting functionality
    "test_vector_vs_scalar.jl",
    "test_multi_emitter.jl",   # Multiple emitters
    "test_spline_psf.jl",      # Spline interpolation
    "test_normalization.jl",   # Normalization checks
]

# Run all test files and collect results
results = Dict{String,Bool}()
for file in test_files
    results[file] = run_test_file(file)
end

# Print summary
println("\n\n")
println("="^50)
println("Test Summary")
println("="^50)

success_count = count(values(results))
total_count = length(results)
println("$success_count/$total_count tests passed")

if success_count == total_count
    println("✓ All tests passed!")
else
    println("✗ Some tests failed:")
    for (file, success) in results
        if !success
            println("  • $file")
        end
    end
end

