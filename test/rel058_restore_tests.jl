@testset "0.5.8 restoration review" verbose = true begin
    M = MicroscopePSFs
    root = dirname(@__DIR__)

    @testset "Review 7a: documented constructor examples" begin
        # Run the actual manual example in its own namespace.
        manual = read(joinpath(root, "docs/src/zernike.md"), String)
        example = split(split(manual, "## Usage in PSF Models")[2], "```julia\n")[2]
        example = first(split(example, "```"))
        sandbox = Module(:ZernikeManualExample)
        Core.eval(sandbox, :(using MicroscopePSFs))
        Base.include_string(sandbox, example)
        @test Core.eval(sandbox, :(psf isa ScalarPSF))
        @test Core.eval(sandbox, :(vector_psf isa VectorPSF))
        @test Core.eval(sandbox, :(psf.zernike_coeffs.phase[4])) == 0.5
        for file in ("dev/test_Scalar3D.jl", "dev/test_auto_diff.jl")
            @test !occursin(r"\bcoeffs\s*=", read(joinpath(root, file), String))
        end
    end

    @testset "Review 7b: Noll example names" begin
        overview = read(joinpath(root, "api_overview.md"), String)
        expected = Dict(
            5 => (2, -2), 6 => (2, 2), 7 => (3, -1),
            11 => (4, 0), 22 => (6, 0)
        )
        names = Dict(
            5 => "Oblique astigmatism", 6 => "Vertical astigmatism",
            7 => "Vertical coma", 11 => "Primary spherical", 22 => "Secondary spherical"
        )
        for line in split(overview, '\n')
            mode = match(r"zc.phase\[(\d+)\].*# (.*)", line)
            isnothing(mode) && continue
            j = parse(Int, mode[1])
            haskey(names, j) || continue
            # A generic astigmatism comment is also valid for j=6.
            occursin("Add astigmatism", mode[2]) && j == 6 && continue
            @test strip(mode[2]) == names[j]
            @test M.noll2nl(j) == expected[j]
        end
    end

    @testset "Review 7c: supplied vector pupil grid" begin
        for N in (16, 32), requested in (16, 32, 64), rotating in (false, true)
            base = PupilFunction(1.4, 0.6, 1.33, ZernikeCoefficients(1); grid_size = N)
            args = rotating ? () : (DipoleVector(1.0, 0.0, 0.0),)
            psf = VectorPSF(1.4, 0.6, args...; base_pupil = base, grid_size = requested)
            reference = VectorPSF(1.4, 0.6, args...; base_pupil = base, grid_size = N)
            @test psf.base_pupil === base
            @test all(p -> size(p.Ex.field) == (N, N), psf.vector_pupils)
            @test dipole_fields(psf, 0.1, -0.2, 0.3) ==
                dipole_fields(reference, 0.1, -0.2, 0.3)
        end

        # Independent matched-index Debye field: transmission is one and the
        # apodization is 1/sqrt(cos(theta)). Do not use production pupil helpers.
        function aberrated_vector_reference(base, dipoles, xyz)
            N = size(base, 1)
            step = 4 / (N - 1) # NA=1, wavelength=0.5, so kmax=2.
            center = (N + 1) / 2
            fields = [zeros(ComplexF64, 2, N, N) for _ in dipoles]
            phase = zeros(ComplexF64, N, N)
            for i in 1:N, j in 1:N
                kx, ky = step * (i - center), step * (j - center)
                radius² = kx^2 + ky^2
                radius² > 4 && continue
                angle = atan(ky, kx)
                c, s = cos(angle), sin(angle)
                sinθ = sqrt(radius²) / 3 # n/wavelength=3.
                cosθ = sqrt(1 - radius² / 9)
                for (field, (px, py, pz)) in zip(fields, dipoles)
                    field[1, j, i] = base[j, i] / sqrt(cosθ) * (
                        px * (cosθ * c^2 + s^2) + py * (cosθ - 1) * c * s +
                            pz * sinθ * c
                    )
                    field[2, j, i] = base[j, i] / sqrt(cosθ) * (
                        px * (cosθ - 1) * c * s + py * (cosθ * s^2 + c^2) +
                            pz * sinθ * s
                    )
                end
                phase[j, i] = cis(
                    2π * (
                        kx * xyz[1] + ky * xyz[2] +
                            sqrt(9 - radius²) * xyz[3]
                    )
                )
            end
            energy = sum(sum(abs2, field) for field in fields) * step^2
            pupils = [field / sqrt(energy) for field in fields]
            values = [
                [sum(field[c, :, :] .* phase) * step^2 for c in 1:2]
                    for field in pupils
            ]
            return pupils, values
        end

        for N in (16, 32), requested in (16, 32, 64), rotating in (false, true)
            coords = range(-1, 1; length = N)
            field = [
                (0.8 + 0.2x + 0.1y^2) * cis(0.3x^2 + 0.2x * y - 0.17y)
                    for y in coords, x in coords
            ]
            base = PupilFunction(1.0, 0.5, 1.5, field)
            args = rotating ? () : (DipoleVector(1.0, 2.0, 3.0),)
            dipoles = rotating ? [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)] :
                [Tuple([1, 2, 3] / sqrt(14))]
            psf = VectorPSF(
                1.0, 0.5, args...; n_medium = 1.5, n_coverslip = 1.5,
                n_immersion = 1.5, base_pupil = base, grid_size = requested
            )
            for xyz in ((0.0, 0.0, 0.0), (0.1, -0.2, 0.3))
                pupils, values = aberrated_vector_reference(field, dipoles, xyz)
                for (actual, expected) in zip(dipole_fields(psf, xyz...), values)
                    @test actual ≈ expected rtol = 1.0e-12 atol = 1.0e-13
                end
                @test psf(xyz...) ≈ sum(sum(abs2, value) for value in values) rtol = 1.0e-12
                if iszero(xyz[3])
                    for (pupil, expected) in zip(psf.vector_pupils, pupils)
                        @test pupil.Ex.field ≈ expected[1, :, :] rtol = 1.0e-12
                        @test pupil.Ey.field ≈ expected[2, :, :] rtol = 1.0e-12
                    end
                end
            end
        end
    end

    @testset "Review 7d: ScalarPSF constructor documentation" begin
        text = string(@doc ScalarPSF(::Real, ::Real, ::Real))
        @test occursin("Create a scalar 3D PSF", text)
        @test occursin("grid_size", text)
        @test occursin("psf(x, y, z)", text)
    end

    @testset "Review 7e: zero positional pupil" begin
        zero_pupil = PupilFunction(1.4, 0.6, 1.518, zeros(ComplexF64, 16, 16))
        for args in ((), (nothing,), (ZernikeCoefficients(1),))
            @test_throws ArgumentError ScalarPSF(1.4, 0.6, 1.518, zero_pupil, args...)
        end
        @test_throws r"scale=0.0.*> 0" ScalarPSF(1.4, 0.6, 1.518, zero_pupil)
        @test_throws ArgumentError ScalarPSF(1.4, 0.6, 1.518; pupil = zero_pupil)
        pupil = PupilFunction(1.4, 0.6, 1.518, ones(ComplexF64, 16, 16))
        psf = ScalarPSF(1.4, 0.6, 1.518, pupil)
        @test psf.pupil === pupil
        @test all(isone, pupil.field)
    end

    @testset "Note 8a: coherent multi-emitter field documentation" begin
        for text in (
                string(
                    @doc integrate_pixels_amplitude(
                        ::AbstractPSF, ::AbstractCamera, ::Vector{<:AbstractEmitter}
                    )
                ),
                string(
                    @doc integrate_pixels_amplitude(
                        ::AbstractPSF, ::AbstractVector, ::AbstractVector,
                        ::Vector{<:AbstractEmitter}
                    )
                ),
            )
            @test occursin("coherent sum", text)
            @test occursin("mutually incoherent", text)
            @test occursin("abs2", text)
            @test occursin("integrate_pixels", text)
        end
        psf = GaussianPSF(0.13)
        camera = IdealCamera(5, 5, 0.1)
        emitter = Emitter2D(0.25, 0.25, 1.0)
        field = integrate_pixels_amplitude(psf, camera, emitter)
        fields = integrate_pixels_amplitude(psf, camera, [emitter, emitter])
        @test fields == 2field
        @test abs2.(fields) ≈ 4abs2.(field)
        @test integrate_pixels(psf, camera, [emitter, emitter]) ==
            2integrate_pixels(psf, camera, emitter)
    end

    @testset "Note 8b: displayed pupil orientation" begin
        manual = read(joinpath(root, "docs/src/zernike.md"), String)
        @test occursin("yreversed=true", manual)
        @test occursin("positive y points down", manual)
        @test occursin("increases clockwise", manual)
        @test !occursin("Posiive y axis is up", manual)
    end

    # Read source examples so the tests also verify the documented calls themselves.
    function documented_block(path, marker)
        text = split(read(joinpath(root, path), String), marker; limit = 2)[2]
        return match(r"```julia\r?\n(.*?)```"s, text)[1]
    end
    function run_documented(code)
        sandbox = Module(gensym(:DocExample))
        Core.eval(sandbox, :(using MicroscopePSFs))
        Base.include_string(sandbox, code)
        return sandbox
    end

    @testset "Follow-up 1: vector aberration example" begin
        example = documented_block("docs/src/psfs/vectorpsf.md", "## Examples")
        sandbox = run_documented(example)
        @test Core.eval(sandbox, :(psf_aberrated isa VectorPSF))
        @test Core.eval(sandbox, :(psf_aberrated.zernike_coeffs.phase[11])) == 0.5
        @test Core.eval(sandbox, :(isfinite(psf_aberrated(0.1, 0.0, 0.0))))
    end

    @testset "Follow-up 2: positional scalar aperture bound" begin
        for T in (Float32, Float64), (na, n) in ((1.4, 1.33), (1.6, 1.52))
            pupil = PupilFunction(T(na), T(0.6), T(n), ones(Complex{T}, 8, 8))
            original = copy(pupil.field)
            for args in ((), (nothing,), (ZernikeCoefficients(1; T),))
                error = try
                    ScalarPSF(T(na), T(0.6), T(n), pupil, args...)
                catch caught
                    caught
                end
                @test error isa ArgumentError
                if error isa ArgumentError
                    @test occursin("nₐ=$(T(na))", error.msg)
                    @test occursin("n=$(T(n))", error.msg)
                    @test occursin("<=", error.msg)
                end
            end
            @test pupil.field == original
            @test_throws ArgumentError ScalarPSF(T(na), T(0.6), T(n); pupil)
        end
        for na in (1.0, 1.33)
            pupil = PupilFunction(na, 0.6, 1.33, ones(ComplexF64, 8, 8))
            @test ScalarPSF(na, 0.6, 1.33, pupil).pupil === pupil
            @test isfinite(ScalarPSF(na, 0.6, 1.33, pupil)(0.0, 0.0, 0.0))
        end
        changelog = read(joinpath(root, "CHANGELOG.md"), String)
        unreleased = first(split(split(changelog, "## [Unreleased]")[2], "## [0.5.7]"))
        @test occursin("NA above the medium index", unreleased)
    end

    @testset "Follow-up 3: scalar docstring aberration example" begin
        example = documented_block("src/psfs/scalar3d.jl", "# Examples")
        sandbox = run_documented(example)
        @test Core.eval(sandbox, :(psf isa ScalarPSF))
        @test Core.eval(sandbox, :(psf.zernike_coeffs.phase[5])) == 0.5
        @test Core.eval(sandbox, :(isfinite(psf(0.1, 0.0, 0.0))))
    end

    @testset "Follow-up 4: documented Noll term count" begin
        overview = read(joinpath(root, "api_overview.md"), String)
        @test occursin("# Create Zernike coefficients (15 Noll terms)", overview)
        @test !occursin("max radial order 15", overview)
        example = documented_block("api_overview.md", "### PSF with Aberrations")
        sandbox = run_documented(example)
        @test Core.eval(sandbox, :(length(zc))) == 15
        @test Core.eval(sandbox, :(psf_aberrated isa ScalarPSF))
    end

    @testset "Doc audit: conventions coefficient setter" begin
        text = read(joinpath(root, "docs/src/conventions.md"), String)
        expression = match(r"`(zc\[[^`]*= 0.5)`", text)[1]
        sandbox = run_documented("zc = ZernikeCoefficients(15)\n" * expression)
        @test Core.eval(sandbox, :(zc.phase[5])) == 0.5
        @test occursin("The Noll indexing scheme is used by default", text)
        @test !occursin("The OSA indexing scheme is used by default", text)
    end

    @testset "Doc audit: experimental HDF5 spline example" begin
        example = documented_block(
            "docs/src/psfs/splinepsf.md",
            "## Working with Experimental PSFs"
        )
        mktempdir() do dir
            path = joinpath(dir, "measured.h5")
            # A non-square field detects reversed x/y dimension assignments.
            values = [
                exp(-i^2 / 5 - j^2 / 8 - k^2 / 6)
                    for i in -2:2, j in -3:3, k in -2:2
            ]
            h5open(path, "w") do file
                file["psf"] = values
            end
            example = replace(example, "path/to/measured_psf.h5" => path)
            sandbox = Module(gensym(:MeasuredPSFExample))
            Core.eval(sandbox, :(using MicroscopePSFs))
            # HDF5 is a package dependency, accessed here without another test dependency.
            Core.eval(sandbox, :(using MicroscopePSFs.HDF5))
            example = replace(example, "using HDF5" => "using MicroscopePSFs.HDF5")
            Base.include_string(sandbox, example)
            spline = Core.eval(sandbox, :(experimental_psf))
            @test spline isa SplinePSF
            @test length(spline.x_range) == 7
            @test length(spline.y_range) == 5
            @test spline.original_grid == values
            @test spline(0.0, 0.0, 0.0) ≈ 1
        end
    end

    @testset "Doc audit: Airy conversion signature" begin
        text = read(joinpath(root, "src/psfs/airy2d.jl"), String)
        signature = match(r"    AiryPSF\(psf::GaussianPSF[^\n]*", text).match
        expression = replace(
            strip(signature),
            "psf::GaussianPSF" => "GaussianPSF(0.15)",
            "λ::Real=0.532" => "λ=0.532", ", λ)" => ", 0.532)"
        )
        sandbox = run_documented("psf = " * expression)
        @test Core.eval(sandbox, :(psf isa AiryPSF))
        @test Core.eval(sandbox, :(psf.λ)) == 0.532
    end

    @testset "Doc audit: tetrapod coefficient example" begin
        example = documented_block("api_overview.md", "### Example 2:")
        sandbox = run_documented(first(split(example, "# Evaluate PSF")))
        @test Core.eval(sandbox, :(length(zc))) >= 22
        @test Core.eval(sandbox, :(psf_tetrapod.zernike_coeffs.phase[22])) == 0.8
        @test Core.eval(sandbox, :(isfinite(psf_tetrapod(0.1, 0.0, 0.0))))
    end

    @testset "Doc audit: width analysis example" begin
        example = documented_block("docs/src/interface.md", "### Using the Same Code")
        sandbox = run_documented(example)
        @test Core.eval(sandbox, :(results_gaussian)) ≈ 2sqrt(2log(2)) * 0.15 atol = 0.03
        @test Core.eval(sandbox, :(0 < results_airy < results_gaussian))
    end

    @testset "Doc audit: astigmatism coefficient example" begin
        example = documented_block("api_overview.md", "### 2. 3D")
        sandbox = run_documented(first(split(example, "# Emitter")))
        @test Core.eval(sandbox, :(length(zc))) >= 6
        @test Core.eval(sandbox, :(psf.zernike_coeffs.phase[6])) == 1
    end

    @testset "Doc audit: broadcast PSF evaluation" begin
        example = documented_block("api_overview.md", "### PSF Evaluation")
        setup = "psf = GaussianPSF(0.15)\nx, y, z = 0.1, 0.2, 0.3\n" *
            "x_vec, y_vec, z_vec = [0.0, 0.1], [0.0, 0.2], [0.0, 0.3]\n"
        sandbox = run_documented(setup * example)
        @test Core.eval(sandbox, :(length(intensity))) == 2
        @test Core.eval(
            sandbox, :(
                intensity ==
                    [psf_3d(x_vec[i], y_vec[i], z_vec[i]) for i in eachindex(x_vec)]
            )
        )
        @test Core.eval(sandbox, :(isfinite(abs2(amp))))
    end

    @testset "Doc audit: OSA conversion example result" begin
        text = read(joinpath(root, "src/zernike/indexing.jl"), String)
        example = match(r"julia> (nl2osa\(4, 2\))\n(\d+)", text)
        sandbox = Module(gensym(:OSAExample))
        Core.eval(sandbox, :(using MicroscopePSFs.Zernike))
        result = Core.eval(sandbox, Meta.parse(example[1]))
        @test result == parse(Int, example[2])
        @test result == 13
    end

    @testset "Doc audit: Gaussian constructor signature" begin
        example = documented_block(
            "docs/src/psfs/gaussianpsf.md",
            "## Constructor and Parameters"
        )
        type_name = match(r"GaussianPSF\(σ::(\w+)\)", example)[1]
        argument_type = Core.eval(Main, Meta.parse(type_name))
        @test hasmethod(GaussianPSF, Tuple{argument_type})
        @test GaussianPSF(0.15) isa GaussianPSF
    end
end

@testset "Review 15: exact coordinates and Gaussian range" verbose = true begin
    @testset "Review 15.1: exact primal window membership" begin
        M = MicroscopePSFs
        @test Base.get_extension(M, :MicroscopePSFsForwardDiffExt) !== nothing

        function window_models(T, width)
            optics = (
                n_medium = T(1.33), n_coverslip = T(1.52),
                n_immersion = T(1.52), grid_size = 17,
            )
            wavelength = width / T(4)
            return (
                ScalarPSF(T(1), wavelength, T(1.33); grid_size = 17),
                VectorPSF(
                    T(1), wavelength, DipoleVector(T(0.3), T(0.4), T(0.5)); optics...
                ),
                VectorPSF(T(1), wavelength; optics...),
            )
        end

        function check_window(psf, x, axis, halfwidth)
            # High-precision distances give an independent membership reference.
            outside = abs(BigFloat(x)) >= BigFloat(halfwidth)
            xyz = ntuple(i -> i == axis ? x : i == 3 ? 0.5 : 0.1, 3)
            value = psf(xyz...)
            @test isfinite(value)
            @test iszero(value) == outside
            @test all(field -> all(iszero, field), dipole_fields(psf, xyz...)) == outside
            if psf isa ScalarPSF || length(psf.vector_pupils) == 1
                @test all(iszero, amplitude(psf, xyz...)) == outside
            end
            return outside, xyz, value
        end

        setprecision(BigFloat, 256) do
            for T in (Float32, Float64), constant in (π, ℯ)
                target = T(constant)
                for width in (prevfloat(target), target, nextfloat(target))
                    for psf in window_models(T, width)
                        pupil = psf isa ScalarPSF ? psf.pupil : first(psf.vector_pupils)
                        halfwidth = M._alias_free_halfwidth(pupil)
                        coordinates = (
                            constant, target, prevfloat(target), nextfloat(target),
                            BigFloat(constant), Rational{BigInt}(BigFloat(constant)),
                        )
                        for axis in 1:2, x in coordinates
                            check_window(psf, x, axis, halfwidth)
                        end
                    end
                end
            end

            for T in (Float32, Float64), psf in window_models(T, T(2))
                pupil = psf isa ScalarPSF ? psf.pupil : first(psf.vector_pupils)
                halfwidth = M._alias_free_halfwidth(pupil)
                @test halfwidth == 2
                coordinates = (
                    19999999999999999 // 10000000000000000, 2 // 1,
                    20000000000000001 // 10000000000000000,
                    prevfloat(BigFloat(2)), BigFloat(2), nextfloat(BigFloat(2)),
                )
                for axis in 1:2, side in (-1, 1), position in coordinates
                    x = side * position
                    outside, xyz, reference = check_window(psf, x, axis, halfwidth)
                    for seed in (-one(x), one(x))
                        dual = ForwardDiff.Dual{Nothing}(x, seed)
                        nested = ForwardDiff.Dual{Nothing}(dual, one(dual) * seed)
                        for d in (dual, nested)
                            dual_xyz = ntuple(i -> i == axis ? d : xyz[i], 3)
                            value = psf(dual_xyz...)
                            @test iszero(value) == outside
                            @test all(
                                field -> all(iszero, field), dipole_fields(psf, dual_xyz...)
                            ) == outside
                            if !outside
                                primal = ForwardDiff.value(ForwardDiff.value(value))
                                @test primal ≈ reference rtol = 64eps(Float64)
                            end
                        end
                    end
                end
            end
        end

        # The exact review reproductions, without converting either coordinate.
        scalar = ScalarPSF(1.0f0, Float32(π) / 4, 1.5f0; grid_size = 17)
        @test scalar(π, 0.0f0, 0.0f0) > 0
        scalar = ScalarPSF(1.2, 0.6, 1.33; grid_size = 17)
        x = 19999999999999999 // 10000000000000000
        for seed in (-1 // 1, 1 // 1)
            dual = ForwardDiff.Dual{Nothing}(x, seed)
            @test ForwardDiff.value(scalar(dual, 0.1, 0.5)) ≈ scalar(x, 0.1, 0.5)
        end
    end

    @testset "Review 15.2: Float32 Gaussian final rounding" begin
        setprecision(BigFloat, 256) do
            widths = (
                nextfloat(0.0f0), floatmin(Float32), 1.0f-25, 1.0f-20, 0.15f0,
                1.0f19, 1.0f20, 1.0f25, floatmax(Float32),
            )
            for width in widths, ratio in (0.0f0, 1.0f0, 2.0f0, 8.0f0, 10.0f0, 16.0f0)
                psf = GaussianPSF(width)
                x, y = ratio * width, zero(width)
                σ = BigFloat(width)
                reference = exp(-(BigFloat(x)^2 + BigFloat(y)^2) / (2σ^2)) /
                    (2big(π) * σ^2)
                expected = Float32(reference)
                value = @inferred psf(x, y)
                @test value isa Float32
                @test isinf(value) == isinf(expected)
                @test iszero(value) == iszero(expected)
                @test value ≈ expected rtol = 2eps(Float32) atol = eps(0.0f0)
                field = @inferred amplitude(psf, x, y)
                expected_field = Float32(sqrt(reference))
                @test field isa Float32
                @test isinf(field) == isinf(expected_field)
                @test iszero(field) == iszero(expected_field)
                @test field ≈ expected_field rtol = 2eps(Float32) atol = eps(0.0f0)
            end
        end
        @test GaussianPSF(1.0f-20)(1.0f-19, 0.0f0) ≈ 3.069701f17
    end
end

@testset "Review 16: Float32 rounding after final operations" verbose = true begin
    function gaussian_field_reference(psf, x, y)
        σ = BigFloat(psf.σ)
        return exp(-(BigFloat(x)^2 + BigFloat(y)^2) / (4σ^2)) / (sqrt(2big(π)) * σ)
    end
    function airy_field_reference(psf, x, y)
        ν = iszero(psf.ν) || !isfinite(psf.ν) ?
            2big(π) * BigFloat(psf.nₐ) / BigFloat(psf.λ) : BigFloat(psf.ν)
        w² = (BigFloat(x)^2 + BigFloat(y)^2) * ν^2
        # Independent convergent Bessel series at the tested radii (w < 9).
        term = one(w²)
        series = term
        for m in 1:160
            term *= -w² / (4m * (m + 1))
            series += term
        end
        return ν / sqrt(4big(π)) * series
    end
    function scalar_field_reference(psf, x, y, z)
        p = psf isa PupilFunction ? psf : psf.pupil
        N = size(p.field, 1)
        step = 2 * (p.nₐ / p.λ) / (N - 1)
        center = (N + 1) / 2
        limit = (p.nₐ / p.λ)^2
        result = zero(Complex{BigFloat})
        for i in 1:N, j in 1:N
            kx, ky = Float64(step) * (i - center), Float64(step) * (j - center)
            # Exact aperture geometry: kx^2 + ky^2 < kmax^2 in whole numbers.
            (2i - N - 1)^2 + (2j - N - 1)^2 < (N - 1)^2 || continue
            kz = sqrt(complex(BigFloat((p.n / p.λ)^2) - BigFloat(kx)^2 - BigFloat(ky)^2))
            phase = 2big(π) * (BigFloat(x) * kx + BigFloat(y) * ky + BigFloat(z) * kz)
            result += Complex{BigFloat}(p.field[j, i]) * exp(im * phase)
        end
        return result * BigFloat(step)^2
    end
    function raw_scalar(field, wavelength = 0.5f0)
        pupil = PupilFunction(1.0f0, wavelength, 1.5f0, fill(ComplexF32(field), 17, 17))
        return ScalarPSF(1.0f0, wavelength, 1.5f0, pupil)
    end
    setprecision(BigFloat, 256) do
        @testset "Gaussian and Airy amplitudes" begin
            for width in (0.15f0, 1.0f-20, 1.0f-40, 1.0f20, 1.0f25, floatmax(Float32))
                psf = GaussianPSF(width)
                field = @inferred amplitude(psf, 0.0f0, 0.0f0)
                expected = Float32(gaussian_field_reference(psf, 0, 0))
                @test field isa Float32
                @test field ≈ expected rtol = 2eps(Float32) atol = eps(0.0f0)
            end
            for wavelength in (0.5f0, 1.0f-20, 1.0f20)
                psf = AiryPSF(1.0f0, wavelength)
                for x in (0.0f0, wavelength)
                    field = @inferred amplitude(psf, x, 0.0f0)
                    reference = airy_field_reference(psf, x, 0)
                    @test field isa Float32
                    @test field ≈ Float32(reference) rtol = 32eps(Float32) atol = eps(0.0f0)
                    @test (@inferred psf(x, 0.0f0)) ≈ Float32(abs2(reference)) rtol =
                        32eps(Float32) atol = eps(0.0f0)
                end
            end
        end
        @testset "Pupil amplitude final quadrature weight" begin
            for (field, wavelength) in ((1.0f0, 0.5f0), (1.0f37, 1.0f10), (1.0f37, 1.0f21))
                psf = raw_scalar(field, wavelength)
                reference = scalar_field_reference(psf, 0, 0, 0)
                actual = amplitude(psf.pupil, 0.0f0, 0.0f0, 0.0f0)
                @test actual ≈ ComplexF32(reference) rtol = 32eps(Float32)
            end
            for N in (16, 17)
                pupil = PupilFunction(
                    1.5f0, 0.5f0, 1.0f0, fill(ComplexF32(1.0f-38), N, N)
                )
                reference = scalar_field_reference(pupil, 0, 0, -10)
                actual = @inferred amplitude(pupil, 0.0f0, 0.0f0, -10.0f0)
                @test actual isa ComplexF32
                @test isfinite(actual)
                @test actual ≈ ComplexF32(reference) rtol = sqrt(eps(Float32))
            end
        end
    end
end

@testset "Review 14: primal field-window boundaries" verbose = true begin
    M = MicroscopePSFs
    @testset "Integer-valued dual coordinates" begin
        pupil = PupilFunction(1.0, 2.0^60, 1.5, ones(ComplexF64, 17, 17))
        psf = ScalarPSF(pupil.nₐ, pupil.λ, pupil.n, pupil)
        halfwidth = Int(M._alias_free_halfwidth(pupil))
        for axis in 1:2, side in (-1, 1)
            # An Int one unit inside the boundary rounds to it in Float64.
            inside = side * (halfwidth - 1)
            xyz = ntuple(i -> i == axis ? inside : 0, 3)
            @test !iszero(psf(xyz...))
            for seed in (-1, 1)
                dual = ForwardDiff.Dual{Nothing}(inside, seed)
                dual_xyz = ntuple(i -> i == axis ? dual : xyz[i], 3)
                @test ForwardDiff.value(psf(dual_xyz...)) ≈ psf(xyz...)
                value = amplitude(pupil, dual_xyz...)
                primal = complex(
                    ForwardDiff.value(real(value)), ForwardDiff.value(imag(value))
                )
                @test primal ≈ amplitude(pupil, xyz...)
            end
        end
    end
    for T in (Float32, Float64), N in (16, 17)
        scalar = ScalarPSF(T(1.2), T(0.6), T(1.33); grid_size = N)
        optics = (
            n_medium = T(1.33), n_coverslip = T(1.52),
            n_immersion = T(1.52), grid_size = N,
        )
        fixed = VectorPSF(
            T(1.2), T(0.6), DipoleVector(T(1), T(0), T(0)); optics...
        )
        rotating = VectorPSF(T(1.2), T(0.6); optics...)
        for (name, psf) in (("scalar", scalar), ("fixed", fixed), ("rotating", rotating))
            @testset "$T grid $N $name" begin
                pupil = psf === scalar ? scalar.pupil : first(psf.vector_pupils)
                halfwidth = M._alias_free_halfwidth(pupil)
                for axis in 1:2, side in (-one(T), one(T))
                    boundary = side * halfwidth
                    xyz = ntuple(i -> i == axis ? boundary : T(i == 3 ? 0.5 : 0.1), 3)
                    @test iszero(psf(xyz...))
                    for seed in (-one(T), one(T))
                        dual = ForwardDiff.Dual{Nothing}(boundary, seed)
                        dual_xyz = ntuple(i -> i == axis ? dual : xyz[i], 3)
                        @test iszero(psf(dual_xyz...))
                        @test all(
                            field -> all(iszero, field), dipole_fields(psf, dual_xyz...)
                        )
                        if psf !== rotating
                            @test all(iszero, amplitude(psf, dual_xyz...))
                        end
                        if psf === scalar
                            @test iszero(amplitude(pupil, dual_xyz...))
                        end

                        # Nested AD must also classify the innermost primal boundary.
                        nested = ForwardDiff.Dual{Nothing}(dual, one(dual) * seed)
                        nested_xyz = ntuple(i -> i == axis ? nested : xyz[i], 3)
                        @test iszero(psf(nested_xyz...))

                        # The nearest representable interior point remains inside;
                        # neither seed direction may expand or shrink the window.
                        inside = side < 0 ? nextfloat(boundary) : prevfloat(boundary)
                        outside = side < 0 ? prevfloat(boundary) : nextfloat(boundary)
                        for (position, is_outside) in ((inside, false), (outside, true))
                            offset = ForwardDiff.Dual{Nothing}(position, seed)
                            offset_xyz = ntuple(i -> i == axis ? offset : xyz[i], 3)
                            @test iszero(psf(offset_xyz...)) == is_outside
                        end
                    end
                end
            end
        end
    end
end

@testset "Review 13: Float32 range and precision documentation" verbose = true begin
    function scalar_range_reference(pupil, xyz)
        N = size(pupil.field, 1)
        step = 2 * (pupil.nₐ / pupil.λ) / (N - 1)
        center = (N + 1) / 2
        result = 0.0 + 0.0im
        for i in 1:N, j in 1:N
            kx, ky = step * (i - center), step * (j - center)
            radius² = kx^2 + ky^2
            (2i - N - 1)^2 + (2j - N - 1)^2 < (N - 1)^2 || continue    # exact aperture
            kz = sqrt(complex(Float64((pupil.n / pupil.λ)^2) - radius²))
            phase = 2π * (xyz[1] * kx + xyz[2] * ky + xyz[3] * kz)
            result += ComplexF64(pupil.field[j, i]) * exp(im * phase)
        end
        return result * Float64(step)^2
    end

    @testset "Review 13.1: Float32 pupil accumulation range" begin
        for (N, wavelength, value) in (
                (16, 1.0f10, ComplexF32(1.0f37)),
                (17, 1.0f10, ComplexF32(1.0f37)),
                (16, 1.0f10, ComplexF32(3.0f38, 3.0f38)),
                (128, 1.0f21, ComplexF32(1.0f37)),
            )
            pupil = PupilFunction(1.0f0, wavelength, 1.5f0, fill(value, N, N))
            psf = ScalarPSF(1.0f0, wavelength, 1.5f0, pupil)
            points = wavelength > 1.0f20 ? ((0.0f0, 0.0f0, 0.0f0),) :
                ((0.0f0, 0.0f0, 0.0f0), (0.1f0, -0.07f0, 0.12f0) .* wavelength)
            for xyz in points
                reference = scalar_range_reference(pupil, Float64.(xyz))
                field = @inferred amplitude(psf, xyz...)
                value = @inferred psf(xyz...)
                @test field isa ComplexF32
                @test value isa Float32
                @test isfinite(field)
                @test isfinite(value)
                @test field ≈ reference rtol = sqrt(eps(Float32))
                @test value ≈ abs2(reference) rtol = sqrt(eps(Float32))
            end
        end
        pupil = PupilFunction(1.0f0, 1.0f10, 1.5f0, fill(ComplexF32(1.0f37), 16, 16))
        psf = ScalarPSF(1.0f0, 1.0f10, 1.5f0, pupil)
        @test amplitude(psf, 0.0f0, 0.0f0, 0.0f0) ≈ 3.057777829354438e17 rtol =
            8eps(Float32)
        @test psf(0.0f0, 0.0f0, 0.0f0) ≈ 9.350005253691537e34 rtol = 8eps(Float32)
        # Raw pupils also support evanescent samples; retain the complex square root.
        evanescent = PupilFunction(
            1.0f0, 1.0f10, 0.8f0, fill(ComplexF32(1.0f37), 16, 16)
        )
        xyz = (1.0f9, -0.7f9, 1.2f9)
        @test amplitude(evanescent, xyz...) ≈
            scalar_range_reference(evanescent, Float64.(xyz)) rtol = sqrt(eps(Float32))
    end

    @testset "Review 13.2: Float32 Gaussian normalization range" begin
        for width in (1.0f19, 1.0f20), ratio in (0.0f0, 1.0f0, 2.0f0)
            psf = GaussianPSF(width)
            x, y = ratio * width, 0.0f0
            reference = exp(-(Float64(x) / Float64(width))^2 / 2) /
                (2π * Float64(width)^2)
            value = @inferred psf(x, y)
            @test value isa Float32
            @test isfinite(value)
            @test value > 0
            @test value ≈ Float32(reference) rtol = 8eps(Float32) atol = eps(0.0f0)
            @test (@inferred amplitude(psf, x, y)) isa Float32
        end
        @test GaussianPSF(1.0f19)(0.0f0, 0.0f0) ≈ 1.59155f-39 atol = eps(0.0f0)
    end

    @testset "Review 13.3: documented scalar precision rule" begin
        docstring = string(@doc ScalarPSF(::Real, ::Real, ::Real))
        guide = read(joinpath(dirname(@__DIR__), "docs/src/psfs/scalarpsf.md"), String)
        for text in (docstring, guide)
            @test occursin("ScalarPSF{Float32}", text)
            @test occursin("Float64 coefficients or pupil components", text)
        end
        for T in (Float32, Float64)
            pupil = PupilFunction(1.4f0, 0.6f0, 1.518f0, ones(Complex{T}, 16, 16))
            options = (
                (; zernike_coeffs = ZernikeCoefficients(15; T)),
                (; pupil), (; pupil_data = copy(pupil.field)),
            )
            for keywords in options
                psf = ScalarPSF(1.4f0, 0.6f0, 1.518f0; grid_size = 16, keywords...)
                @test psf isa ScalarPSF{T}
                @test (@inferred psf(0.0f0, 0.0f0, 0.0f0)) isa T
                @test (@inferred amplitude(psf, 0.0f0, 0.0f0, 0.0f0)) isa Complex{T}
                @test psf(0, 0.0f0, 0.0f0) isa Float64
            end
        end
    end
end

@testset "Review 12: signed integer field-window boundaries" verbose = true begin
    scalar = ScalarPSF(1.2, 0.6, 1.33; grid_size = 16)
    fixed = VectorPSF(1.2, 0.6, DipoleVector(1.0, 0.0, 0.0); grid_size = 16)
    rotating = VectorPSF(1.2, 0.6; grid_size = 16)
    plane = PlanePSF(scalar; dz = 0.1)
    for T in (Int8, Int16, Int32, Int64)
        @testset "$T" begin
            for extreme in (typemin(T), typemax(T)), axis in 1:2
                xyz = axis == 1 ? (extreme, zero(T), zero(T)) :
                    (zero(T), extreme, zero(T))
                for psf in (scalar, fixed, rotating, plane)
                    @test iszero(psf(xyz...))
                    @test all(field -> all(iszero, field), dipole_fields(psf, xyz...))
                    if psf !== rotating
                        @test all(iszero, amplitude(psf, xyz...))
                    end
                    # A NaN coordinate takes precedence over the lateral window.
                    nan_xyz = (xyz[1], xyz[2], NaN)
                    @test isnan(psf(nan_xyz...))
                    @test all(field -> all(isnan, field), dipole_fields(psf, nan_xyz...))
                    if psf !== rotating
                        @test all(isnan, amplitude(psf, nan_xyz...))
                    end
                end
                @test iszero(amplitude(scalar.pupil, xyz...))
                @test isnan(amplitude(scalar.pupil, xyz[1], xyz[2], NaN))
            end
        end
    end
end

@testset "Review 11: default scalar constructor precision" verbose = true begin
    # Values measured at ef969db with the default grid and no supplied pupil/coefficients.
    cases = (
        (
            "Float16", (Float16(1.4), Float16(0.6), Float16(1.518)),
            17.076213057723137, 0.16854272117576785 + 2.5585121393234274im,
        ),
        (
            "Int/Float32 aperture", (1, 0.5f0, 1.5f0),
            12.54287308574597, -1.1175561829939553 + 2.3088947956471757im,
        ),
        (
            "Int/Float32 wavelength", (1.0f0, 1, 1.5f0),
            3.1357182714364926, 0.8881106210689277 + 1.3768048883784205im,
        ),
        (
            "Int/Float32 index", (1.0f0, 0.5f0, 2),
            12.54287308574597, -2.4672568901659884 + 0.7612718868196726im,
        ),
        (
            "Float64/Float32 aperture", (1.4, 0.6f0, 1.518f0),
            17.072242565475804, 0.16548617485477884 + 2.5589890633420924im,
        ),
        (
            "Float64/Float32 wavelength", (1.4f0, 0.6, 1.518f0),
            17.072243340796028, 0.16548596397054666 + 2.558989085463303im,
        ),
        (
            "Float64/Float32 index", (1.4f0, 0.6f0, 1.518),
            17.07224198400845, 0.16548616987329187 + 2.558989069930216im,
        ),
        (
            "Float32", (1.4f0, 0.6f0, 1.518f0),
            17.07224198400845, 0.16548614163186381 + 2.5589890739426595im,
        ),
        (
            "Float64", (1.4, 0.6, 1.518),
            17.072243922263382, 0.16548602543489196 + 2.5589890708503im,
        ),
    )

    function constructor_energy64(pupil)
        N = size(pupil.field, 1)
        step = 2 * (pupil.nₐ / pupil.λ) / (N - 1)
        center = (N + 1) / 2
        energy = 0.0
        for i in 1:N, j in 1:N
            kx, ky = step * (i - center), step * (j - center)
            if kx^2 + ky^2 < (pupil.nₐ / pupil.λ)^2
                energy += abs2(ComplexF64(pupil.field[j, i]))
            end
        end
        return energy * Float64(step)^2
    end

    for (name, optics, center, field) in cases
        @testset "$name" begin
            all_float32 = all(x -> x isa Float32, optics)
            T = all_float32 ? Float32 : Float64
            tolerance = all_float32 ? sqrt(eps(Float32)) : 32eps(Float64)
            psf = @inferred ScalarPSF(optics...)
            @test psf isa ScalarPSF{T}
            @test psf.pupil isa PupilFunction{T}
            @test size(psf.pupil.field) == (128, 128)
            @test isnothing(psf.zernike_coeffs)
            @test constructor_energy64(psf.pupil) ≈ 1 rtol = max(tolerance, 1.0e-12)
            @test (@inferred psf(0.0, 0.0, 0.0)) ≈ center rtol = tolerance
            @test (@inferred amplitude(psf, 0.1, -0.07, 0.12)) ≈ field rtol = tolerance
            @test (@inferred psf(0.1, -0.07, 0.12)) ≈ abs2(field) rtol = tolerance
            @test (@inferred psf(0.0f0, 0.0f0, 0.0f0)) isa T
            @test (@inferred amplitude(psf, 0.0f0, 0.0f0, 0.0f0)) isa Complex{T}
            @test psf(0.0f0, 0.0f0, 0.0f0) ≈ center rtol = tolerance
            @test amplitude(psf, 0.0f0, 0.0f0, 0.0f0) ≈ sqrt(center) rtol = tolerance
        end
    end
end

@testset "0.5.8 review 10" verbose = true begin
    M = MicroscopePSFs

    # Independent quadrature with the grid arithmetic used at ef969db.
    function reference_samples(pupil)
        N = size(pupil.field, 1)
        step = 2 * (pupil.nₐ / pupil.λ) / (N - 1)
        center = (N + 1) / 2
        samples = [
            (i, j, step * (i - center), step * (j - center))
                for i in 1:N for j in 1:N
                if (2i - N - 1)^2 + (2j - N - 1)^2 < (N - 1)^2    # exact aperture geometry
        ]
        return step, samples
    end
    function reference_amplitude(pupil, x, y, z)
        step, samples = reference_samples(pupil)
        result = zero(ComplexF64)
        for (i, j, kx, ky) in samples
            kz = sqrt(complex((pupil.n / pupil.λ)^2 - kx^2 - ky^2))
            result += pupil.field[j, i] * exp(2π * im * (x * kx + y * ky + z * kz))
        end
        return result * step^2
    end
    function reference_energy(pupil)
        step, samples = reference_samples(pupil)
        return sum(abs2(ComplexF64(pupil.field[j, i])) for (i, j, _, _) in samples) *
            step^2
    end

    @testset "Review 10.1: odd-grid aperture and Float32 quadrature" begin
        tolerance = sqrt(eps(Float32))
        for N in (9, 17, 25, 33)
            psf = ScalarPSF(1.4f0, 0.6f0, 1.518f0; grid_size = N)
            @test reference_energy(psf.pupil) ≈ 1 rtol = tolerance
            if N == 17
                # The four samples exactly on the rim are outside, in every precision.
                @test length(reference_samples(psf.pupil)[2]) == 193
                @test psf(0.0f0, 0.0f0, 0.0f0) ≈ 193 * Float64(M.kpixelsize(psf.pupil))^2 rtol =
                    tolerance
            end
            for xyz in ((0.0f0, 0.0f0, 0.0f0), (0.1f0, -0.07f0, 0.12f0))
                reference = reference_amplitude(psf.pupil, Float64.(xyz)...)
                @test amplitude(psf, xyz...) ≈ reference rtol = tolerance
                @test psf(xyz...) ≈ abs2(reference) rtol = tolerance
                @test amplitude(psf, xyz...) ≈ amplitude(psf, Float64.(xyz)...) rtol =
                    tolerance
                z = ForwardDiff.Dual{Nothing}(xyz[3], 1.0f0)
                dual_amplitude = amplitude(psf, xyz[1], xyz[2], z)
                dual_value = complex(
                    ForwardDiff.value(real(dual_amplitude)),
                    ForwardDiff.value(imag(dual_amplitude))
                )
                @test dual_value ≈ reference rtol = 32eps(Float64)
                @test ForwardDiff.derivative(t -> psf(xyz[1], xyz[2], t), xyz[3]) ≈
                    ForwardDiff.derivative(
                    t -> abs2(reference_amplitude(psf.pupil, xyz[1], xyz[2], t)),
                    xyz[3]
                ) rtol = 64eps(Float64) atol = 64eps(Float64)
            end
            # The aperture and normalization also agree for an arbitrary complex field.
            field = ComplexF32[
                complex(1 + i / N, j / N) for j in 1:N, i in 1:N
            ]
            pupil = PupilFunction(1.4f0, 0.6f0, 1.518f0, field)
            @test M._pupil_energy(pupil) ≈ reference_energy(pupil) rtol = tolerance
            M.normalize!(pupil)
            @test reference_energy(pupil) ≈ 1 rtol = tolerance
            @test amplitude(pupil, 0.0f0, 0.0f0, 0.0f0) ≈
                reference_amplitude(pupil, 0.0, 0.0, 0.0) rtol = tolerance
            # Exactly representable optics allow a separate Float64 PSF comparison.
            single = ScalarPSF(1.0f0, 0.5f0, 1.5f0; grid_size = N)
            double = ScalarPSF(1.0, 0.5, 1.5; grid_size = N)
            @test single(0.1f0, -0.07f0, 0.12f0) ≈
                double(Float64(0.1f0), Float64(-0.07f0), Float64(0.12f0)) rtol =
                tolerance
        end
    end

    @testset "Review 10.2: NaN propagation outside the window" begin
        for T in (Float32, Float64)
            scalar = ScalarPSF(T(1.4), T(0.6), T(1.518); grid_size = 16)
            fixed = VectorPSF(
                T(1.4), T(0.6), DipoleVector(T(1), T(0), T(0)); grid_size = 16
            )
            rotating = VectorPSF(T(1.4), T(0.6); grid_size = 16)
            for nan_axis in 1:3, outside_axis in 1:2, outside in (-20, 20, -Inf, Inf)
                nan_axis == outside_axis && continue
                xyz = ntuple(3) do axis
                    axis == nan_axis ? T(NaN) : axis == outside_axis ? T(outside) : zero(T)
                end
                for psf in (scalar, fixed, rotating)
                    @test isnan(psf(xyz...))
                    if psf !== rotating
                        @test all(isnan, amplitude(psf, xyz...))
                    end
                    @test all(field -> all(isnan, field), dipole_fields(psf, xyz...))
                end
                @test isnan(amplitude(scalar.pupil, xyz...))
            end
            for psf in (scalar, fixed, rotating), outside in (-20, 20, -Inf, Inf)
                @test isnan(psf(T(outside), T(-outside), T(NaN)))
                if psf !== rotating
                    @test all(isnan, amplitude(psf, T(outside), T(-outside), T(NaN)))
                end
            end
        end
    end

    # For inputs outside all-Float32, the ef969db formulas evaluated on Float64 copies of the
    # stored parameter and the coordinates: narrow numbers are widened before squaring.
    function reference_airy(psf, x, y)
        ν, x, y = Float64(psf.ν), Float64(x), Float64(y)
        r² = x^2 + y^2
        w² = r² * ν^2
        w² < 1.0e-10 && return ν / sqrt(4π) * (1 - w² / 8)
        w = sqrt(r²) * ν
        return ν / sqrt(4π) * (2 * besselj1(w) / w)
    end
    function reference_gaussian(psf, x, y)
        σ, x, y = Float64(psf.σ), Float64(x), Float64(y)
        return 1 / (2π * σ^2) * exp(-(x^2 + y^2) / (2σ^2))
    end

    @testset "Review 10.3: Float32 scope and inferred mixed inputs" begin
        for T in (Float16, Float32, Float64)
            airy = AiryPSF(T(1.4), T(0.6))
            gaussian = GaussianPSF(T(0.13))
            pupil = PupilFunction(T(1.4), T(0.6), T(1.518), ones(Complex{T}, 17, 17))
            scalar = ScalarPSF(T(1.4), T(0.6), T(1.518), pupil)
            for X in (Float16, Float32, Float64, Int), Y in (Float32, Float64, Int)
                R = T === Float32 && X === Float32 && Y === Float32 ? Float32 : Float64
                for (x, y) in ((zero(X), zero(Y)), (one(X), zero(Y)))
                    for psf in (airy, gaussian)
                        value = @inferred psf(x, y)
                        field = @inferred amplitude(psf, x, y)
                        @test value isa R
                        @test field isa R
                        if R === Float64
                            reference = if psf isa AiryPSF
                                reference_airy(psf, x, y)
                            else
                                sqrt(reference_gaussian(psf, x, y))
                            end
                            @test field === reference
                            @test value === (
                                psf isa AiryPSF ? abs2(reference) :
                                    reference_gaussian(psf, x, y)
                            )
                        end
                    end
                    xyz = (x, y, zero(Y))
                    value = @inferred scalar(xyz...)
                    field = @inferred amplitude(scalar, xyz...)
                    @test value isa R
                    @test field isa Complex{R}
                    if R === Float64
                        reference = reference_amplitude(pupil, xyz...)
                        @test field ≈ reference rtol = 32eps(Float64)
                        @test value ≈ abs2(reference) rtol = 64eps(Float64)
                    end
                end
                # Early returns have the same inferred type as interior evaluation.
                @test (@inferred amplitude(scalar, X(20), zero(Y), zero(Y))) isa Complex{R}
                @test (@inferred scalar(X(20), zero(Y), zero(Y))) isa R
            end
        end
    end

    @testset "Review 10.4: nonzero positional pupil below energy range" begin
        pupil = PupilFunction(
            1.4f0, 0.6f0, 1.518f0, fill(ComplexF32(1.0f-23), 16, 16)
        )
        @test M._pupil_energy(pupil) == 0.0f0
        original = copy(pupil.field)
        for args in ((), (nothing,), (ZernikeCoefficients(1; T = Float32),))
            psf = ScalarPSF(1.4f0, 0.6f0, 1.518f0, pupil, args...)
            @test psf.pupil === pupil
            @test pupil.field == original
            @test (@inferred psf(0.0, 0.0, 0.0)) ≈ 2.7715255298306554e-44
            @test (@inferred psf(0.0f0, 0.0f0, 0.0f0)) > 0.0f0
        end
        # Nonzero samples outside the aperture do not make a zero pupil usable.
        outside = zeros(ComplexF32, 16, 16)
        outside[1, 1] = 1.0f0
        @test_throws ArgumentError ScalarPSF(
            1.4f0, 0.6f0, 1.518f0, PupilFunction(1.4f0, 0.6f0, 1.518f0, outside)
        )
        for value in (ComplexF32(NaN), ComplexF32(Inf))
            field = copy(original)
            field[8, 8] = value
            @test_throws ArgumentError ScalarPSF(
                1.4f0, 0.6f0, 1.518f0, PupilFunction(1.4f0, 0.6f0, 1.518f0, field)
            )
        end
    end

    @testset "Narrow inputs widened before squaring" begin
        # Before: each square was taken in the narrow type, so these gave NaN, zero, the
        # peak, or a wrapped integer square. References use the stored values exactly.
        setprecision(BigFloat, 256) do
            gaussian_big(σ, x, y) = (s = BigFloat(σ);
                exp(-(BigFloat(x)^2 + BigFloat(y)^2) / (2s^2)) / (2big(π) * s^2))
            function airy_big(ν, x, y)
                v = BigFloat(ν)
                w = sqrt(BigFloat(x)^2 + BigFloat(y)^2) * v
                a = iszero(w) ? v / sqrt(4big(π)) : v / sqrt(4big(π)) * 2besselj1(w) / w
                return abs2(a)
            end
            gaussian_cases = (
                (1.0f-25, 1.0e-24, 0.0),          # Float32 width: square underflowed, NaN
                (1.0f-25, 0.0, 0.0),
                (Float16(1.0e-4), 1.0e-4, 0.0),   # Float16 width: NaN
                (Float16(1.0e-4), 0.0, 0.0),
                (1.0f25, 1.0e25, 0.0),            # Float32 width: square overflowed, zero
                (1.0e-25, 1.0f-24, 0.0f0),        # Float32 position: square underflowed, peak
                (5.0, Int8(12), Int8(0)),         # Int8 position: square wrapped
                (Float16(0.13), 0.1, 0.0),        # Float16 width: square rounded, 2e-4 off
            )
            for (σ, x, y) in gaussian_cases
                psf = GaussianPSF(σ)
                reference = gaussian_big(σ, x, y)
                value = psf(x, y)
                @test value isa Float64
                @test value ≈ reference rtol = 1.0e-13
                @test amplitude(psf, x, y) ≈ sqrt(reference) rtol = 1.0e-13
            end
            airy_cases = (
                (Float16(1.4), Float16(0.02), 0.0),   # Float16 ν² overflowed: NaN at center
                (Float16(1.4), Float16(0.02), 0.01),
                (1.4f0, 1.0f-19, 0.0),                # Float32 ν² overflowed: NaN at center
                (1.4f0, 1.0f-19, 1.0e-20),
                (1.4, 0.6, Int8(12)),                 # Int8 position: square wrapped
            )
            for (na, λ, x) in airy_cases
                psf = AiryPSF(na, λ)
                reference = airy_big(psf.ν, x, zero(x))
                value = psf(x, zero(x))
                @test value isa Float64
                @test value ≈ reference rtol = 1.0e-13
            end
        end
        # Float64 parameters with Float64 positions: the ef969db arithmetic, bit for bit.
        for (x, y) in ((0.0, 0.0), (0.0731, -0.0412), (0.2, 0.35), (1.7, -2.3))
            gaussian = GaussianPSF(0.13)
            @test gaussian(x, y) === 1 / (2π * 0.13^2) * exp(-(x^2 + y^2) / (2 * 0.13^2))
            airy = AiryPSF(1.4, 0.6)
            r² = x^2 + y^2
            w² = r² * airy.ν^2
            expected = w² < 1.0e-10 ? airy.ν / sqrt(4π) * (1 - w² / 8) :
                airy.ν / sqrt(4π) * (2 * besselj1(sqrt(r²) * airy.ν) / (sqrt(r²) * airy.ν))
            @test amplitude(airy, x, y) === expected
        end
    end

    @testset "Integer coordinates into a Float64 image" begin
        # A preallocated Float64 image shows the wrapped values themselves.
        function into_float64(psf, ex, ey, emitter; sampling)
            image = zeros(length(ey) - 1, length(ex) - 1)
            M.integrate_pixels!(image, psf, ex, ey, emitter; sampling, threaded = false)
            return image
        end
        for (psf, ex, ey, x₀, sampling) in (
                (AiryPSF(1.4, 10.0), Int8[0, 12], Int8[0, 12], Int8(0), 1),
                (GaussianPSF(100.0), Int8[100, 120], Int8[-10, 10], Int8(-120), 2),
            )
            emitter = Emitter2D(x₀, zero(x₀), one(x₀))
            reference = into_float64(
                psf, Float64.(ex), Float64.(ey), Emitter2D(Float64(x₀), 0.0, 1.0); sampling
            )
            @test into_float64(psf, ex, ey, emitter; sampling) == reference
            @test all(>(0), reference)
        end
    end

    @testset "Integer pixel edges and positions" begin
        # Before: Int8 arithmetic wrapped. A pixel from 0 to 12 microns had area 144, which
        # is -112 in Int8; an edge at 100 minus an emitter at -120 gave 220, which is -36.
        # Now each matches the same numbers given as Float64, bit for bit.
        f64(v) = Float64.(v)
        cases = (
            (AiryPSF(1.4, 10.0), Int8[0, 12], Int8[0, 12], (Int8(0), Int8(0))),
            (AiryPSF(1.4, 0.6), Int8[-3, -1, 1, 3], Int8[-3, -1, 1, 3], (Int8(1), Int8(-1))),
            (GaussianPSF(100.0), Int8[100, 120], Int8[-10, 10], (Int8(-120), Int8(0))),
            (GaussianPSF(0.5), Int16[-2, 0, 2], Int16[-2, 0, 2], (Int16(0), Int16(1))),
            (ScalarPSF(1.2, 0.6, 1.33), Int8[-2, 0, 2], Int8[-2, 0, 2], (Int8(0), Int8(0))),
        )
        for (psf, ex, ey, (x₀, y₀)) in cases, sampling in (1, 2, 3)
            # An emitter holds one number type, so integer emitters count integer photons.
            I = typeof(x₀)
            emitter = psf isa ScalarPSF ? Emitter3D(x₀, y₀, zero(I), I(100)) :
                Emitter2D(x₀, y₀, I(100))
            reference_emitter = psf isa ScalarPSF ?
                Emitter3D(Float64(x₀), Float64(y₀), 0.0, 100.0) :
                Emitter2D(Float64(x₀), Float64(y₀), 100.0)
            image = integrate_pixels(psf, ex, ey, emitter; sampling, threaded = false)
            reference = integrate_pixels(
                psf, f64(ex), f64(ey), reference_emitter; sampling, threaded = false
            )
            @test image == reference
            @test eltype(image) == Float64
            @test all(>=(0), image)
            pair = integrate_pixels(psf, ex, ey, [emitter, emitter]; sampling, threaded = false)
            @test pair == integrate_pixels(
                psf, f64(ex), f64(ey), [reference_emitter, reference_emitter];
                sampling, threaded = false
            )
            field = integrate_pixels_amplitude(psf, ex, ey, emitter; sampling, threaded = false)
            @test field == integrate_pixels_amplitude(
                psf, f64(ex), f64(ey), reference_emitter; sampling, threaded = false
            )
        end
    end

    @testset "Integer support radius and optical parameters" begin
        # Before: an Int8 emitter at 100 plus an Int8 radius of 40 wrapped to -116, so the
        # support region was empty and the image zero; AiryPSF(1, 1) and GaussianPSF(1)
        # failed with method errors.
        ex, ey = Int8[95, 100, 105], Int8[-5, 0, 5]
        emitter = Emitter2D(Int8(100), Int8(0), Int8(100))
        reference_emitter = Emitter2D(100.0, 0.0, 100.0)
        for psf in (GaussianPSF(0.5), AiryPSF(1.4, 0.6)), support in (Int8(40), Int8(3))
            image = integrate_pixels(psf, ex, ey, emitter; support, threaded = false)
            reference = integrate_pixels(
                psf, Float64.(ex), Float64.(ey), reference_emitter;
                support = Float64(support), threaded = false
            )
            @test image == reference
            @test sum(reference) > 0
            field = integrate_pixels_amplitude(psf, ex, ey, emitter; support, threaded = false)
            @test field == integrate_pixels_amplitude(
                psf, Float64.(ex), Float64.(ey), reference_emitter;
                support = Float64(support), threaded = false
            )
            pair = integrate_pixels(psf, ex, ey, [emitter, emitter]; support, threaded = false)
            @test pair == integrate_pixels(
                psf, Float64.(ex), Float64.(ey), [reference_emitter, reference_emitter];
                support = Float64(support), threaded = false
            )
        end
        airy = AiryPSF(1, 1)
        @test airy isa AiryPSF{Float64}
        @test airy(0.3, 0.1) === AiryPSF(1.0, 1.0)(0.3, 0.1)
        gaussian = GaussianPSF(1)
        @test gaussian isa GaussianPSF{Float64}
        @test gaussian(0.3, 0.1) === GaussianPSF(1.0)(0.3, 0.1)
        @test ScalarPSF(1, 1, 2) isa ScalarPSF{Float64}
        dipole = DipoleVector(0.0, 0.0, 1.0)
        @test VectorPSF(1, 1, dipole; n_medium = 2, n_coverslip = 2, n_immersion = 2) isa
            VectorPSF{Float64}
    end

    @testset "Integer axial coordinates" begin
        # Before: the 4Pi PSF negated an unsigned z, -UInt8(7) == 249; integration passed an
        # integer z through unconverted, so a Float32 PlanePSF computed 7 + dz in Float32.
        s = ScalarPSF(1.2, 0.6, 1.4; grid_size = 9)
        fourpi = FourPiPSF(s, s; phase = 0.3)
        s32 = ScalarPSF(1.2f0, 0.6f0, 1.4f0; grid_size = 9)
        plane = PlanePSF(s32; dz = 0.02f0)
        ex, ey = UInt8[0, 1, 5], UInt8[0, 2, 5]
        emitter = Emitter3D(UInt8(0), UInt8(0), UInt8(7), UInt8(100))
        reference_emitter = Emitter3D(0.0, 0.0, 7.0, 100.0)
        for psf in (fourpi, plane), sampling in (1, 2)
            image = integrate_pixels(psf, ex, ey, emitter; sampling, threaded = false)
            @test image == integrate_pixels(
                psf, Float64.(ex), Float64.(ey), reference_emitter; sampling, threaded = false
            )
            pair = integrate_pixels(psf, ex, ey, [emitter, emitter]; sampling, threaded = false)
            @test pair == integrate_pixels(
                psf, Float64.(ex), Float64.(ey), [reference_emitter, reference_emitter];
                sampling, threaded = false
            )
            @test psf(0.5, 0.5, UInt8(7)) == psf(0.5, 0.5, 7.0)
        end
        mc = MultiChannelPSF(; channels = [s32, s32], z_offset = [0.0f0, 0.2f0])
        @test channel_position(mc, 2, UInt8(3), UInt8(4), UInt8(7)) ==
            channel_position(mc, 2, 3.0, 4.0, 7.0)
    end

    @testset "Spline intensity cut off at zero" begin
        # Before: a spline of positive samples interpolated to -0.28 at x = 0.396 and gave
        # negative photon counts there.
        r = -1.0:0.25:1.0
        gaussian = GaussianPSF(0.13)
        spline = SplinePSF(gaussian, r, r)
        @test spline(0.396, 0.0) == 0.0
        grid = -1.0:0.004:1.0
        @test all(>=(0), spline(x, y) for x in grid, y in (0.0, 0.1, 0.396))
        image = integrate_pixels(
            spline, [0.395, 0.397], [-0.001, 0.001], Emitter2D(0.0, 0.0, 1.0);
            sampling = 1, threaded = false
        )
        @test all(>=(0), image)
        @test real(amplitude(spline, 0.396, 0.0)) == 0 && imag(amplitude(spline, 0.396, 0.0)) == 0
        # Samples are still reproduced, and NaN still propagates.
        for x in (-0.5, 0.0, 0.25), y in (-0.25, 0.0)
            @test spline(x, y) ≈ gaussian(x, y) rtol = 1.0e-12 atol = 1.0e-14
        end
        @test isnan(spline(NaN, 0.0))
        scalar = ScalarPSF(1.2, 0.6, 1.33; grid_size = 17)
        spline3 = SplinePSF(scalar, r, r, -0.5:0.25:0.5)
        @test all(>=(0), spline3(x, y, z) for x in -1.0:0.01:1.0, y in (0.0, 0.3), z in (-0.4, 0.1))
        @test isnan(spline3(0.0, 0.0, NaN))
    end

    @testset "Integer photons with single-precision positions" begin
        # Before: the Gaussian route multiplied an integer photon count by Float32 pixel
        # fractions, so the image held Float32 products, unlike photons given as Float64.
        @eval Main mutable struct MixedPhotonEmitter{P} <: MicroscopePSFs.AbstractEmitter
            x::Float32
            y::Float32
            photons::P
        end
        ex, ey = Float32[-0.1, 0.05, 0.25], Float32[-0.15, 0.0, 0.2]
        for psf in (GaussianPSF(0.13f0), AiryPSF(1.4f0, 0.6f0)), I in (Int, Int8, UInt16)
            whole = Main.MixedPhotonEmitter(0.01f0, -0.02f0, I(100))
            floating = Main.MixedPhotonEmitter(0.01f0, -0.02f0, 100.0)
            image = integrate_pixels(psf, ex, ey, whole; threaded = false)
            @test image == integrate_pixels(psf, ex, ey, floating; threaded = false)
            @test eltype(image) == Float64
        end
    end

    @testset "Integer counts into a Float32 image; all-integer vector optics" begin
        # Before: a preallocated Float32 image was scaled by an integer count in Float32, so
        # Int32(16777217) rounded before the multiplication; VectorPSF with every optical
        # argument an integer, z_stage included, failed with a TypeError.
        ex, ey = Float32[-0.1, 0.05, 0.25], Float32[-0.15, 0.0, 0.2]
        for psf in (AiryPSF(1.4f0, 0.6f0), GaussianPSF(0.13f0))
            into(photons) = (image = zeros(Float32, 2, 2);
                M.integrate_pixels!(image, psf, ex, ey,
                    Main.MixedPhotonEmitter(0.01f0, -0.02f0, photons); threaded = false);
                image)
            @test into(Int32(16777217)) == into(16777217.0)
            @test into(Int64(10)^12) == into(1.0e12)
        end
        vector = VectorPSF(1, 1; n_medium = 2, n_coverslip = 2, n_immersion = 2, z_stage = 0,
            grid_size = 9)
        @test vector isa VectorPSF{Float64}
        @test PupilFunction(1, 1, 2, ones(Complex{Int}, 9, 9)) isa PupilFunction{Float64}
    end

    @testset "Integer Zernike coefficients; supplied pupil of another precision" begin
        # Before: integer Zernike coefficients kept an integer element type, and evaluating
        # the pupil threw InexactError; ScalarPSF with Float64 optics and a Float32 pupil
        # threw a MethodError, since no conversion between pupil types existed.
        zi = ZernikeCoefficients([1, 0, 0, 1], [0, 0, 0, 1])
        zf = ZernikeCoefficients([1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0, 1.0])
        @test PupilFunction(1, 1, 2, zi; grid_size = 9).field ==
            PupilFunction(1.0, 1.0, 2.0, zf; grid_size = 9).field
        si = ScalarPSF(1, 1, 2; zernike_coeffs = zi, grid_size = 9)
        @test si isa ScalarPSF{Float64}
        @test si(0.1, 0.2, 0.0) == ScalarPSF(1.0, 1.0, 2.0; zernike_coeffs = zf, grid_size = 9)(0.1, 0.2, 0.0)
        p32 = PupilFunction(1.4f0, 0.6f0, 1.518f0, ones(ComplexF32, 9, 9))
        p64 = PupilFunction(Float64(1.4f0), Float64(0.6f0), Float64(1.518f0), ComplexF64.(p32.field))
        # The positional route keeps a supplied pupil, so the widened Float32 pupil is the
        # Float64 one exactly. The keyword route normalizes a copy in the pupil's own
        # precision, so the two agree to Float32 rounding.
        positional = ScalarPSF(1.4, 0.6, 1.518, p32)
        @test positional isa ScalarPSF{Float64}
        @test positional.pupil.field == p64.field
        @test positional(0.1, 0.2, 0.3) == ScalarPSF(1.4, 0.6, 1.518, p64)(0.1, 0.2, 0.3)
        keyword = ScalarPSF(1.4, 0.6, 1.518; pupil = p32)
        @test keyword isa ScalarPSF{Float64}
        @test keyword(0.1, 0.2, 0.3) ≈ ScalarPSF(1.4, 0.6, 1.518; pupil = p64)(0.1, 0.2, 0.3) rtol = 1.0e-6
        @test convert(PupilFunction{Float32}, p32) === p32
        @test ScalarPSF(big"1.4", big"0.6", big"1.518"; pupil = p64) isa ScalarPSF{BigFloat}
    end

    @testset "Supplied pupil left unchanged; aperture membership exact in every precision" begin
        # Before: the keyword ScalarPSF normalized the caller's pupil in place; on an odd grid
        # a Float32 pupil counted the four samples exactly on the aperture rim and the same
        # pupil in Float64 did not, so their normalizations differed by 4 percent.
        mk(T) = PupilFunction(T(1.4f0), T(0.6f0), T(1.518f0), ones(Complex{T}, 9, 9))
        p = mk(Float64)
        before = copy(p.field)
        ScalarPSF(1.4, 0.6, 1.518; pupil = p)
        @test p.field == before
        inside(q) = count(M._pupil_sample(q, i, j)[4] for i in 1:9, j in 1:9)
        @test inside(mk(Float32)) == inside(mk(Float64)) == 45
        for T in (Float32, Float64)
            vector = VectorPSF(T(1.2), T(0.6), DipoleVector(T(0), T(0), T(1));
                n_medium = T(1.4), n_coverslip = T(1.4), n_immersion = T(1.4), grid_size = 9)
            v = vector.vector_pupils[1]
            c = (9 + 1) / 2
            @test count(M._vector_pupil_sample(i, j, c, M.kpixelsize(v), M.kmax(v)^2)[4]
                for i in 1:9, j in 1:9) == 49
        end
        a = ScalarPSF(1.4, 0.6, 1.518; pupil = mk(Float32))
        b = ScalarPSF(1.4, 0.6, 1.518; pupil = mk(Float64))
        for xyz in ((0.1, 0.2, 0.3), (0.0, 0.0, 0.0), (0.4, -0.3, 0.5))
            @test a(xyz...) ≈ b(xyz...) rtol = 1.0e-5
        end
        vf = [VectorPSF(T(1.2), T(0.6), DipoleVector(T(0.3), T(0.4), T(sqrt(0.75)));
            n_medium = T(1.4), n_coverslip = T(1.4), n_immersion = T(1.4), grid_size = 9)
            for T in (Float32, Float64)]
        for xyz in ((0.1, 0.2, 0.3), (0.0, 0.0, 0.0), (0.4, -0.3, 0.5))
            @test vf[1](xyz...) ≈ vf[2](xyz...) rtol = 1.0e-4
        end
    end
end
