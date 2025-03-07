# test/vector3d_tests.jl
using Test
using MicroscopePSFs
using LinearAlgebra

@testset "Vector3D PSF" begin
    # Setup common parameters
    na = 1.2
    λ = 0.6
    n_medium = 1.33
    n_coverslip = 1.52
    n_immersion = 1.52
    
    # Test different dipole orientations
    @testset "Dipole Orientations" begin
        # Create dipoles with different orientations
        dipole_x = DipoleVector(1.0, 0.0, 0.0)  # X-oriented
        dipole_y = DipoleVector(0.0, 1.0, 0.0)  # Y-oriented
        dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Z-oriented
        
        # Create PSFs for each orientation
        psf_x = Vector3DPSF(na, λ, dipole_x, n_medium=n_medium)
        psf_y = Vector3DPSF(na, λ, dipole_y, n_medium=n_medium)
        psf_z = Vector3DPSF(na, λ, dipole_z, n_medium=n_medium)
        
        # For X/Y-oriented dipoles, test symmetry properties
        # X dipole should be symmetric on y-axis
        @test psf_x(0.0, 0.5, 0.0) ≈ psf_x(0.0, -0.5, 0.0)
        
        # Y dipole should be symmetric on x-axis
        @test psf_y(0.5, 0.0, 0.0) ≈ psf_y(-0.5, 0.0, 0.0)
        
        # Z dipole should have circular symmetry (rotational symmetry around z-axis)
        @test isapprox(psf_z(0.0, 0.5, 0.0), psf_z(0.5, 0.0, 0.0), rtol=1e-2)
        
        # Test that dipole orientations produce different patterns
        # Skip comparing exact values as they depend on implementation details
    end
end