using Logging: with_logger

@testset "0.5.8 sweep" verbose = true begin
    M = MicroscopePSFs

    @testset "Should 1: documented Zernike radians" begin
        for file in ("docs/src/zernike.md", "docs/src/psfs/scalarpsf.md")
            text = read(joinpath(@__DIR__, "..", file), String)
            @test !occursin(r"waves RMS|λ RMS|units of wavelength", text)
            @test occursin("rad RMS", text)
        end
        zc = ZernikeCoefficients(15)
        zc.phase[4] = 0.5
        p = PupilFunction(1.4, 0.6, 1.518, zc; grid_size = 65)
        @test angle(p.field[33, 65]) ≈ 0.5sqrt(3)
    end

    @testset "Should 2: scalar grid size" begin
        @test size(ScalarPSF(1.4, 0.6, 1.518).pupil.field) == (128, 128)
        for grid_size in (32, 64, 128)
            psf = ScalarPSF(1.4, 0.6, 1.518; grid_size = grid_size)
            @test size(psf.pupil.field) == (grid_size, grid_size)
        end
        pupil = PupilFunction(1.4, 0.6, 1.518, ZernikeCoefficients(1); grid_size = 32)
        @test ScalarPSF(1.4, 0.6, 1.518; pupil = pupil, grid_size = 64).pupil === pupil
        @test size(
            ScalarPSF(
                1.4, 0.6, 1.518; pupil_data = ones(ComplexF64, 32, 32),
                grid_size = 64
            ).pupil.field
        ) == (32, 32)
    end

    @testset "Should 3: exact Gaussian pixels on every route" begin
        psf = GaussianPSF(0.13)
        ex = [-0.15, -0.01, 0.1, 0.25, 0.42]
        ey = [-0.2, 0.0, 0.12, 0.3]
        camera = IdealCamera(ex, ey)
        emitters = [Emitter2D(0.1, 0.05, 1000.0), Emitter2D(0.0, 0.12, 230.0)]
        function exact(e)
            return [
                e.photons / 4 *
                    (
                    erf((ex[i + 1] - e.x) / (sqrt(2) * psf.σ)) -
                        erf((ex[i] - e.x) / (sqrt(2) * psf.σ))
                ) *
                    (
                    erf((ey[j + 1] - e.y) / (sqrt(2) * psf.σ)) -
                        erf((ey[j] - e.y) / (sqrt(2) * psf.σ))
                )
                    for j in 1:(length(ey) - 1), i in 1:(length(ex) - 1)
            ]
        end
        expected = sum(exact, emitters)
        for sampling in (1, 2, 3, 8), threaded in (false, true)
            for e in emitters
                reference = exact(e)
                @test integrate_pixels(psf, camera, e; sampling, threaded) ≈
                    reference rtol = 1.0e-12
                @test integrate_pixels(psf, ex, ey, e; sampling, threaded) ≈
                    reference rtol = 1.0e-12
                result = similar(reference)
                M.integrate_pixels!(result, psf, camera, e; sampling, threaded)
                @test result ≈ reference rtol = 1.0e-12
                M.integrate_pixels!(result, psf, ex, ey, e; sampling, threaded)
                @test result ≈ reference rtol = 1.0e-12
                support = (ex[2], ex[4], ey[2], ey[4])
                cropped = zeros(size(reference))
                cropped[2:3, 2:3] .= reference[2:3, 2:3]
                @test integrate_pixels(psf, ex, ey, e; support, sampling, threaded) ≈
                    cropped rtol = 1.0e-12
            end
            @test integrate_pixels(psf, camera, emitters; sampling, threaded) ≈
                expected rtol = 1.0e-12
            @test integrate_pixels(psf, ex, ey, emitters; sampling, threaded) ≈
                expected rtol = 1.0e-12
        end
    end

    @testset "Should 5: spline z warning once" begin
        xr = -0.5:0.1:0.5
        sp3 = SplinePSF([exp(-x^2 - y^2 - z^2) for y in xr, x in xr, z in xr], xr, xr, xr)
        camera = IdealCamera(20, 20, 0.1)
        e = Emitter3D(1.0, 1.0, 0.8, 1000.0)
        logger = Test.TestLogger(; respect_maxlog = true)
        with_logger(logger) do
            @test all(iszero, integrate_pixels(sp3, camera, e))
            @test all(
                iszero, integrate_pixels(
                    sp3, camera.pixel_edges_x,
                    camera.pixel_edges_y, e
                )
            )
            @test all(iszero, integrate_pixels(sp3, camera, [e, e]))
        end
        @test length(logger.logs) == 1
        if length(logger.logs) == 1
            message = string(only(logger.logs).message)
            @test occursin("0.8", message)
            @test occursin("-0.5", message)
            @test occursin("0.5", message)
        end
        @test_logs integrate_pixels(sp3, camera, Emitter3D(1.0, 1.0, 0.5, 1.0))
    end

    @testset "Should 7: photon-count documentation" begin
        text = string(
            @doc integrate_pixels(
                ::AbstractPSF, ::AbstractCamera,
                ::AbstractEmitter
            )
        )
        @test !occursin("Values are normalized to sum to 1", text)
        @test occursin("Values are photon counts", text)
        image = integrate_pixels(
            GaussianPSF(0.13), IdealCamera(20, 20, 0.1),
            Emitter2D(1.0, 1.0, 1234.0)
        )
        @test sum(image) ≈ 1234.0 atol = 1.0e-9
    end

    @testset "Should 9: inferred Float32 Gaussian" begin
        psf = GaussianPSF(0.13f0)
        for x in (0.0f0, 0.1f0)
            @test (@inferred psf(x, 0.0f0)) isa Float32
            @test (@inferred amplitude(psf, x, 0.0f0)) isa Float32
        end
    end
end
