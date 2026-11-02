using Logging: with_logger
using ForwardDiff

@testset "0.5.8 sweep" verbose = true begin
    M = MicroscopePSFs

    @testset "AiryPSF center curvature" begin
        psf = AiryPSF(1.4, 0.533)
        hessian = ForwardDiff.hessian(v -> psf(v[1], v[2]), zeros(2))
        curvature = -(psf.ν^2 / (4π)) * psf.ν^2 / 2
        @test hessian[1, 1] ≈ curvature rtol = 1.0e-12
        @test hessian[2, 2] ≈ curvature rtol = 1.0e-12
        @test hessian[1, 2] == 0
        @test hessian[2, 1] == 0
    end

    @testset "AiryPSF values outside the center series" begin
        psf = AiryPSF(1.4, 0.533)
        coordinates = range(-0.5, 0.5; length = 11)
        for x in coordinates, y in coordinates
            iszero(x) && iszero(y) && continue
            # Preserve origin/main's arithmetic outside the center series.
            r = sqrt(x^2 + y^2)
            w = r * psf.ν
            reference = psf.ν / sqrt(4π) * (2 * besselj1(w) / w)
            @test amplitude(psf, x, y) === reference
            @test psf(x, y) === abs2(reference)
        end
    end

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

    @testset "Review 1: 2D emitters on 3D splines" begin
        xr = -0.5:0.1:0.5
        spline = SplinePSF(
            [exp(-x^2 - y^2 - z^2) for y in xr, x in xr, z in xr], xr, xr, xr
        )
        edges = [-0.2, 0.0, 0.2]
        camera = IdealCamera(edges, edges)
        emitter = Emitter2D(0.0, 0.0, 1000.0)
        reference = integrate_pixels(
            spline, edges, edges, Emitter3D(0.0, 0.0, 0.0, 1000.0); threaded = false
        )
        @test sum(reference) ≈ 156.065 atol = 0.0005 rtol = 0
        for threaded in (false, true)
            image = @test_logs integrate_pixels(spline, edges, edges, emitter; threaded)
            @test image == reference
            result = similar(reference)
            @test_logs M.integrate_pixels!(result, spline, edges, edges, emitter; threaded)
            @test result == reference
            @test_logs M.integrate_pixels!(result, spline, camera, emitter; threaded)
            @test result == reference
            batch = @test_logs integrate_pixels(
                spline, edges, edges, [emitter, emitter];
                threaded
            )
            @test batch == 2reference
            batch = @test_logs integrate_pixels(
                spline, camera, [emitter, emitter]; threaded
            )
            @test batch == 2reference
            # The allocating single-emitter camera route requires an explicit z.
            @test_throws ArgumentError integrate_pixels(spline, camera, emitter; threaded)
        end
    end

    @testset "Review 1: Float32 Gaussian tails" begin
        psf = GaussianPSF(0.13f0)
        emitter = Emitter2D(0.0f0, 0.0f0, 1000.0f0)
        tail_edges = [0.75f0, 0.85f0]
        center_edges = [-0.05f0, 0.05f0]
        for threaded in (false, true), (ex, ey) in (
                    (tail_edges, center_edges), (center_edges, tail_edges),
                )
            image = integrate_pixels(psf, ex, ey, emitter; threaded)
            @test eltype(image) == Float32
            @test only(image) ≈ 1.18312f-6 rtol = 2.0f-5
            camera = IdealCamera(ex, ey)
            @test integrate_pixels(psf, camera, emitter; threaded) == image
            result = similar(image)
            M.integrate_pixels!(result, psf, ex, ey, emitter; threaded)
            @test result == image
            @test integrate_pixels(psf, ex, ey, [emitter, emitter]; threaded) == 2image
        end
    end
end
