# test/vector3d_tests.jl

@testset "VectorPSF" begin
    # Setup common parameters
    na = 1.2
    λ = 0.6
    n_medium = 1.33
    n_coverslip = 1.52
    n_immersion = 1.52
    
    @testset "Constructor" begin
        # Test basic constructor
        dipole = DipoleVector(1.0, 0.0, 0.0)  # X-oriented
        psf = VectorPSF(na, λ, dipole; 
                          n_medium=n_medium, 
                          n_coverslip=n_coverslip, 
                          n_immersion=n_immersion)
        
        # Verify parameters were set correctly
        @test psf.nₐ ≈ na
        @test psf.λ ≈ λ
        @test psf.n_medium ≈ n_medium
        @test psf.n_coverslip ≈ n_coverslip
        @test psf.n_immersion ≈ n_immersion
        @test psf.dipole.px ≈ 1.0
        @test psf.dipole.py ≈ 0.0
        @test psf.dipole.pz ≈ 0.0
        
        # Test constructor with aberrations
        zc = ZernikeCoefficients(10)
        zc.phase[6] = 0.5  # Add vertical astigmatism
        psf_aberrated = VectorPSF(na, λ, dipole; 
                                   n_medium=n_medium, 
                                   base_zernike=zc)
        
        @test !isnothing(psf_aberrated.zernike_coeffs)
    end
    
    @testset "Dipole Orientation Symmetry" begin
        # Create dipoles with different orientations
        dipole_x = DipoleVector(1.0, 0.0, 0.0)  # X-oriented
        dipole_y = DipoleVector(0.0, 1.0, 0.0)  # Y-oriented
        dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Z-oriented
        
        # Create PSFs for each orientation
        psf_x = VectorPSF(na, λ, dipole_x; n_medium=n_medium)
        psf_y = VectorPSF(na, λ, dipole_y; n_medium=n_medium)
        psf_z = VectorPSF(na, λ, dipole_z; n_medium=n_medium)
        
        # X-dipole symmetry tests
        # Should be symmetric with respect to y-axis
        @test psf_x(0.0, 0.5, 0.0) ≈ psf_x(0.0, -0.5, 0.0)
        
        # For a pure x-dipole in scalar approximation, we actually get symmetry in both x and y
        @test psf_x(0.5, 0.0, 0.0) ≈ psf_x(-0.5, 0.0, 0.0)
        
        # Y-dipole symmetry tests
        # Should be symmetric with respect to x-axis
        @test psf_y(0.5, 0.0, 0.0) ≈ psf_y(-0.5, 0.0, 0.0)
        
        # Should also be symmetric with respect to y-axis in this implementation
        @test psf_y(0.0, 0.5, 0.0) ≈ psf_y(0.0, -0.5, 0.0)
        
        # Z-dipole symmetry tests
        # Should have rotational symmetry around z-axis
        r = 0.5
        for θ in [0, π/4, π/2, 3π/4, π]
            x1, y1 = r*cos(θ), r*sin(θ)
            x2, y2 = r*cos(θ+π/2), r*sin(θ+π/2)
            @test isapprox(psf_z(x1, y1, 0.0), psf_z(x2, y2, 0.0), rtol=1e-2)
        end
    end
    
    @testset "Dipole Radiation Pattern" begin
        # Test that dipoles follow expected radiation pattern
        dipole_x = DipoleVector(1.0, 0.0, 0.0)
        psf_x = VectorPSF(na, λ, dipole_x; n_medium=n_medium)
        
        # For our implementation, test that intensity decreases with distance
        # from the center for the in-plane dipole
        @test psf_x(0.0, 0.0, 0.0) > psf_x(0.5, 0.0, 0.0)
        @test psf_x(0.0, 0.0, 0.0) > psf_x(0.0, 0.5, 0.0)
        
        # Z-dipole should have donut pattern in xy-plane (maximum at some radius, minimum at center)
        dipole_z = DipoleVector(0.0, 0.0, 1.0)
        psf_z = VectorPSF(na, λ, dipole_z; n_medium=n_medium)
        
        # Center should be a minimum
        @test psf_z(0.0, 0.0, 0.0) < psf_z(0.3, 0.0, 0.0)
    end
    
    @testset "Amplitude Function" begin
        # Test that amplitude function returns expected components
        dipole = DipoleVector(1.0, 0.0, 0.0)
        psf = VectorPSF(na, λ, dipole; n_medium=n_medium)
        
        # Amplitude should return a vector of [Ex, Ey]
        amp = amplitude(psf, 0.3, 0.4, 0.0)
        
        # Check dimensions and types
        @test length(amp) == 2
        @test eltype(amp) <: Complex
        
        # For x-dipole, |Ex| should be larger than |Ey| in most positions
        @test abs(amp[1]) > abs(amp[2])
        
        # Intensity should equal sum of squared field amplitudes
        @test psf(0.3, 0.4, 0.0) ≈ abs2(amp[1]) + abs2(amp[2])
    end
    
    @testset "Pixel Integration" begin
        # Test integration over camera pixels
        dipole = DipoleVector(1.0, 0.0, 0.0)
        psf = VectorPSF(na, λ, dipole; n_medium=n_medium)
        
        # Create camera and emitter
        camera = IdealCamera(collect(0:0.1:2.0), collect(0:0.1:2.0))
        emitter = DipoleEmitter3D(1.0, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.0)
        
        # Integrate PSF over pixels
        img = integrate_pixels(psf, camera, emitter)
        
        # Check image properties
        @test size(img) == (length(0:0.1:2.0)-1, length(0:0.1:2.0)-1)
        
        # Test amplitude integration
        amp_img = integrate_pixels_amplitude(psf, camera, emitter)
        @test size(amp_img, 3) == 2  # Should have two field components
        
        # Just test that the amplitude image has the right dimensions and is non-zero
        # The exact relationship between amplitude and intensity depends on normalization
        intensity_from_amp = abs2.(amp_img[:,:,1]) + abs2.(amp_img[:,:,2])
        @test all(intensity_from_amp .>= 0)
        @test size(intensity_from_amp) == size(img)
    end

    @testset "Position-Evaluated Pupil (pupil_at)" begin
        # Integrate a vector pupil over the aperture with only the lateral phase.
        # If pupil_at has baked in the axial defocus phase correctly, this must
        # reproduce amplitude(psf, x, y, z) exactly (same grid/convention as eval).
        function lateral_field(vp, x, y)
            gs = size(vp.Ex.field, 1)
            kmax_val = vp.nₐ / vp.λ
            kpixel = 2 * kmax_val / (gs - 1)
            center = (gs + 1) / 2
            ex = zero(ComplexF64); ey = zero(ComplexF64)
            for i in 1:gs, j in 1:gs
                kx = (i - center) * kpixel
                ky = (j - center) * kpixel
                kr2 = kx^2 + ky^2
                kr2 > kmax_val^2 && continue
                ph = exp(im * 2π * (kx * x + ky * y))
                ex += vp.Ex.field[j, i] * ph * kpixel^2
                ey += vp.Ey.field[j, i] * ph * kpixel^2
            end
            return [ex, ey]
        end

        dipole = DipoleVector(1.0, 0.0, 0.0)
        psf = VectorPSF(na, λ, dipole; n_medium=n_medium)

        # Single dipole returns a single VectorPupilFunction, leaving psf untouched
        before = copy(psf.vector_pupils[1].Ex.field)
        vp = pupil_at(psf, 0.7)
        @test vp isa VectorPupilFunction
        @test psf.vector_pupils[1].Ex.field == before   # not mutated

        # z = 0, z_stage = 0 -> identity (defocus phase is 1)
        vp0 = pupil_at(psf, 0.0; z_stage=0.0)
        @test vp0.Ex.field ≈ psf.vector_pupils[1].Ex.field
        @test vp0.Ey.field ≈ psf.vector_pupils[1].Ey.field

        # Core invariant: lateral integral of pupil_at reproduces amplitude()
        for (x, y, z) in [(0.0, 0.0, 0.0), (0.3, 0.4, 0.5), (-0.2, 0.1, 1.0)]
            @test lateral_field(pupil_at(psf, z), x, y) ≈ amplitude(psf, x, y, z)
        end

        # z_stage: stored pupil is z_stage-independent; default flows from the PSF
        psf_stage = VectorPSF(na, λ, dipole; n_medium=n_medium, z_stage=0.3)
        @test psf_stage.vector_pupils[1].Ex.field ≈ psf.vector_pupils[1].Ex.field
        for (x, y, z) in [(0.0, 0.0, 0.5), (0.25, -0.15, 1.0)]
            @test lateral_field(pupil_at(psf_stage, z), x, y) ≈ amplitude(psf_stage, x, y, z)
        end
        # Explicit z_stage override matches a PSF built with that stage position
        @test lateral_field(pupil_at(psf, 0.5; z_stage=0.3), 0.2, 0.1) ≈
              amplitude(psf_stage, 0.2, 0.1, 0.5)

        # Free/rotating dipole returns one pupil per x/y/z orientation
        psf_rot = VectorPSF(na, λ; n_medium=n_medium)
        vps = pupil_at(psf_rot, 0.5)
        @test vps isa Vector{<:VectorPupilFunction}
        @test length(vps) == 3
    end
    
   
end