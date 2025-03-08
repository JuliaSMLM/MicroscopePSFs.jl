# test/io_tests.jl

@testset "IO Functions" begin
    # Create temporary directory for test files
    temp_dir = mktempdir()
    
    @testset "Gaussian2D" begin
        # Create and save a Gaussian2D PSF
        σ = 0.15
        psf = Gaussian2D(σ)
        filename = joinpath(temp_dir, "gaussian2d.h5")
        
        # Test saving with metadata
        metadata = Dict("description" => "Test Gaussian PSF", "sigma_um" => σ)
        save_psf(filename, psf, metadata=metadata)
        
        # Test loading
        loaded_psf = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_psf isa Gaussian2D
        @test loaded_psf.σ ≈ psf.σ
        
        # Verify function evaluation
        for _ in 1:5
            x, y = rand() * 0.3 - 0.15, rand() * 0.3 - 0.15
            @test loaded_psf(x, y) ≈ psf(x, y)
        end
    end
    
    @testset "Airy2D" begin
        # Create and save an Airy2D PSF
        nₐ = 1.4
        λ = 0.532
        psf = Airy2D(nₐ, λ)
        filename = joinpath(temp_dir, "airy2d.h5")
        
        # Save PSF
        save_psf(filename, psf)
        
        # Load PSF
        loaded_psf = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_psf isa Airy2D
        @test loaded_psf.nₐ ≈ psf.nₐ
        @test loaded_psf.λ ≈ psf.λ
        @test loaded_psf.ν ≈ psf.ν
        
        # Verify function evaluation
        for _ in 1:5
            x, y = rand() * 0.3 - 0.15, rand() * 0.3 - 0.15
            @test loaded_psf(x, y) ≈ psf(x, y)
        end
    end
    
    @testset "Scalar3DPSF" begin
        # Create and save a Scalar3DPSF
        nₐ = 1.2
        λ = 0.6
        n = 1.33
        
        # With Zernike coefficients
        zc = ZernikeCoefficients(8)
        # Add some astigmatism
        zc.mag[6] = 0.2  # Z_6 = astigmatism
        
        psf = Scalar3DPSF(nₐ, λ, n, coeffs=zc)
        filename = joinpath(temp_dir, "scalar3d.h5")
        
        # Save PSF
        save_psf(filename, psf)
        
        # Load PSF
        loaded_psf = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_psf isa Scalar3DPSF
        @test loaded_psf.nₐ ≈ psf.nₐ
        @test loaded_psf.λ ≈ psf.λ
        @test loaded_psf.n ≈ psf.n
        
        # Verify Zernike coefficients were saved/loaded correctly
        @test !isnothing(loaded_psf.zernike_coeffs)
        @test loaded_psf.zernike_coeffs.mag[6] ≈ 0.2
        
        # Verify function evaluation
        for _ in 1:5
            x, y, z = rand() * 0.3 - 0.15, rand() * 0.3 - 0.15, rand() * 0.4 - 0.2
            @test loaded_psf(x, y, z) ≈ psf(x, y, z)
        end
    end
    
    @testset "SplinePSF" begin
        # Create a SplinePSF from an analytic PSF
        gauss_psf = Gaussian2D(0.15)
        x_range = y_range = range(-1.0, 1.0, length=21)
        z_range = range(-0.5, 0.5, length=5)
        
        # Create a 3D SplinePSF
        spline_psf = SplinePSF(gauss_psf, x_range, y_range, z_range)
        filename = joinpath(temp_dir, "spline_psf.h5")
        
        # Save PSF
        save_psf(filename, spline_psf)
        
        # Load PSF
        loaded_psf = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_psf isa SplinePSF
        @test length(loaded_psf.x_range) == length(spline_psf.x_range)
        @test first(loaded_psf.x_range) ≈ first(spline_psf.x_range)
        @test last(loaded_psf.x_range) ≈ last(spline_psf.x_range)
        @test length(loaded_psf.z_range) == length(spline_psf.z_range)
        
        # Verify function evaluation
        for _ in 1:5
            x, y, z = rand() * 1.8 - 0.9, rand() * 1.8 - 0.9, rand() * 0.9 - 0.45
            @test isapprox(loaded_psf(x, y, z), spline_psf(x, y, z), rtol=1e-3)
        end
        
        # Create and test a 2D SplinePSF
        spline_psf_2d = SplinePSF(gauss_psf, x_range, y_range)
        filename_2d = joinpath(temp_dir, "spline_psf_2d.h5")
        
        # Save and load
        save_psf(filename_2d, spline_psf_2d)
        loaded_psf_2d = load_psf(filename_2d)
        
        # Verify 2D parameters
        @test loaded_psf_2d.z_range === nothing
        
        # Verify function evaluation for 2D
        for _ in 1:5
            x, y = rand() * 1.8 - 0.9, rand() * 1.8 - 0.9
            @test isapprox(loaded_psf_2d(x, y), spline_psf_2d(x, y), rtol=1e-3)
        end
    end
    
    @testset "Vector3DPSF" begin
        # Create a Vector3DPSF
        nₐ = 1.2
        λ = 0.6
        n_medium = 1.33
        n_coverslip = 1.52
        n_immersion = 1.52
        
        # Create dipole orientation
        dipole = DipoleVector(1.0, 0.0, 0.0)  # x-oriented dipole
        
        # Create with Zernike coefficients
        zc = ZernikeCoefficients(10)
        # Add some coma
        zc.mag[8] = 0.15  # Z_8 = coma
        
        psf = Vector3DPSF(nₐ, λ, dipole; 
                         n_medium=n_medium,
                         n_coverslip=n_coverslip,
                         n_immersion=n_immersion,
                         coeffs=zc)
        
        filename = joinpath(temp_dir, "vector3d.h5")
        
        # Save PSF
        save_psf(filename, psf)
        
        # Load PSF
        loaded_psf = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_psf isa Vector3DPSF
        @test loaded_psf.nₐ ≈ psf.nₐ
        @test loaded_psf.λ ≈ psf.λ
        @test loaded_psf.n_medium ≈ psf.n_medium
        @test loaded_psf.n_coverslip ≈ psf.n_coverslip
        @test loaded_psf.n_immersion ≈ psf.n_immersion
        
        # Verify dipole orientation
        @test loaded_psf.dipole.px ≈ psf.dipole.px
        @test loaded_psf.dipole.py ≈ psf.dipole.py
        @test loaded_psf.dipole.pz ≈ psf.dipole.pz
        
        # Verify Zernike coefficients were saved/loaded correctly
        @test !isnothing(loaded_psf.zernike_coeffs)
        @test loaded_psf.zernike_coeffs.mag[8] ≈ 0.15
        
        # Verify function evaluation
        for _ in 1:5
            x, y, z = rand() * 0.3 - 0.15, rand() * 0.3 - 0.15, rand() * 0.4 - 0.2
            @test isapprox(loaded_psf(x, y, z), psf(x, y, z), rtol=1e-3)
        end
    end
    
    @testset "PupilFunction" begin
        # Create a pupil function
        nₐ = 1.2
        λ = 0.6
        n = 1.33
        
        # Create a simple field
        grid_size = 64
        field = zeros(ComplexF64, grid_size, grid_size)
        for i in 1:grid_size, j in 1:grid_size
            r = sqrt((i - grid_size/2)^2 + (j - grid_size/2)^2) / (grid_size/2)
            if r < 1.0
                field[i, j] = exp(-2 * r^2) * exp(1im * r^2)
            end
        end
        
        pupil = PupilFunction(nₐ, λ, n, field)
        filename = joinpath(temp_dir, "pupil.h5")
        
        # Save pupil
        save_psf(filename, pupil)
        
        # Load pupil
        loaded_pupil = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_pupil isa PupilFunction
        @test loaded_pupil.nₐ ≈ pupil.nₐ
        @test loaded_pupil.λ ≈ pupil.λ
        @test loaded_pupil.n ≈ pupil.n
        
        # Verify field values
        @test size(loaded_pupil.field) == size(pupil.field)
        @test all(isapprox.(loaded_pupil.field, pupil.field))
    end
    
    @testset "ZernikeCoefficients" begin
        # Create Zernike coefficients
        zc = ZernikeCoefficients(15)
        # Set some values
        zc.mag[4] = 0.1  # Defocus
        zc.mag[5] = 0.2  # Astigmatism
        zc.phase[6] = 0.5
        
        filename = joinpath(temp_dir, "zernike.h5")
        
        # Save Zernike coefficients
        save_psf(filename, zc)
        
        # Load Zernike coefficients
        loaded_zc = load_psf(filename)
        
        # Verify type and parameters
        @test loaded_zc isa ZernikeCoefficients
        @test length(loaded_zc.mag) == length(zc.mag)
        
        # Verify coefficient values
        @test loaded_zc.mag[4] ≈ zc.mag[4]
        @test loaded_zc.mag[5] ≈ zc.mag[5]
        @test loaded_zc.phase[6] ≈ zc.phase[6]
    end
    
    # Clean up temporary files
    rm(temp_dir, recursive=true)
end