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
                if (step * (i - center))^2 + (step * (j - center))^2 <
                (pupil.nₐ / pupil.λ)^2
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
                @test length(reference_samples(psf.pupil)[2]) == 197
                @test psf(0.0f0, 0.0f0, 0.0f0) ≈ 16.75867889171337 rtol = tolerance
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

    # These formulas preserve the ef969db arithmetic for inputs outside all-Float32.
    function reference_airy(psf, x, y)
        r² = x^2 + y^2
        w² = r² * psf.ν^2
        w² < 1.0e-10 && return psf.ν / sqrt(4π) * (1 - w² / 8)
        w = sqrt(r²) * psf.ν
        return psf.ν / sqrt(4π) * (2 * besselj1(w) / w)
    end
    function reference_gaussian(psf, x, y)
        return 1 / (2π * psf.σ^2) * exp(-(x^2 + y^2) / (2psf.σ^2))
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
end
