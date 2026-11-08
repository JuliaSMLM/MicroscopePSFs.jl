using Test, MicroscopePSFs, HDF5

@testset "0.5.7 review regressions" verbose=true begin
    M = MicroscopePSFs
    energy(vp) = (sum(abs2, vp.Ex.field) + sum(abs2, vp.Ey.field)) *
                 M.kpixelsize(vp)^2
    vector(args...; kwargs...) = VectorPSF(1.2, 0.68, args...;
        n_medium=1.33, n_coverslip=1.33, n_immersion=1.33,
        grid_size=16, kwargs...)
    points = [(0.0, 0.0, 0.0), (0.1, -0.2, 0.3)]

    @testset "I1 joint update normalization" begin
        zc = ZernikeCoefficients(15)
        zc.phase[6] = 0.3
        for args in [(), (DipoleVector(1.0, 0.0, 0.0),)]
            p = vector(args...; base_zernike=zc)
            before = [p(x...) for x in points]
            fields = [(copy(v.Ex.field), copy(v.Ey.field)) for v in p.vector_pupils]
            M.update_pupils!(p)
            @test sum(energy, p.vector_pupils) ≈ 1
            @test [p(x...) for x in points] ≈ before
            for (v, (ex, ey)) in zip(p.vector_pupils, fields)
                @test v.Ex.field ≈ ex
                @test v.Ey.field ≈ ey
            end
            zc.phase[6] += 0.1
            M.update_pupils!(p)
            fresh = vector(args...; base_zernike=zc)
            @test [p(x...) for x in points] ≈ [fresh(x...) for x in points]
        end
    end

    @testset "I2 coefficient sum defines the pupil" begin
        zc = ZernikeCoefficients(15)
        zc.mag[6] = 0.2
        zc.phase[4] = 0.3
        @test PupilFunction(1.2, 0.68, 1.33, zc; grid_size=17).field ==
              M.Zernike.evaluate_pupil(zc, 17)
        flat = PupilFunction(1.2, 0.68, 1.33, ZernikeCoefficients(1); grid_size=17)
        @test flat.field[9, 9] == 1
        @test flat.field[1, 1] == 0
    end

    @testset "I3 radial polynomials through order 40" begin
        # Independent exact factorial sum, evaluated with 256-bit arithmetic.
        oracle(n, m, r) = sum((-1)^k * factorial(big(n-k)) /
            (factorial(big(k)) * factorial(big((n+m)÷2-k)) *
             factorial(big((n-m)÷2-k))) * BigFloat(r)^(n-2k)
            for k in 0:(n-m)÷2)
        for n in 0:40, m in n:-2:0
            @test radialpolynomial(n, m, 1.0) ≈ 1 atol=2e-13
            for r in (0.0, 0.125, 0.5, 0.9, 0.99)
                @test radialpolynomial(n, m, r) ≈ Float64(oracle(n, m, r)) atol=3e-13
            end
        end
        @test radialpolynomial(2, 0, 0.5) == -0.5
        @test radialpolynomial(4, 0, 0.5) == -0.125
        @test radialpolynomial(3, 1, 0.5) == -0.625
        @test radialpolynomial(40, 0, 1.1) == 0
    end

    @testset "I4 pupil geometry and aperture energy" begin
        @test_throws ArgumentError PupilFunction(1.2, 0.68, 1.33, ones(ComplexF64, 7, 5))
        @test_throws ArgumentError PupilFunction(1.2, 0.68, 1.33, ones(ComplexF64, 1, 1))
        for bad in (0.0, Inf, NaN)
            p = PupilFunction(1.2, 0.68, 1.33, fill(complex(bad), 8, 8))
            @test_throws ArgumentError M.normalize!(p)
        end
        p = PupilFunction(1.2, 0.68, 1.33, ones(ComplexF64, 32, 32))
        kpix = M.kpixelsize(p)
        mask = [(2i - 33)^2 + (2j - 33)^2 < 31^2    # exact aperture of the 32-sample grid
                for j in 1:32, i in 1:32]
        M.normalize!(p)
        @test sum(abs2, p.field[mask]) * kpix^2 ≈ 1
        p.field[.!mask] .= 1e6
        before = copy(p.field[mask])
        M.normalize!(p)
        @test p.field[mask] ≈ before
        outside = PupilFunction(1.2, 0.68, 1.33, complex.(.!mask))
        @test_throws ArgumentError M.normalize!(outside)
    end

    @testset "I5 grazing angle samples" begin
        for args in [(), (DipoleVector(1.0, 0.0, 0.0),),
                     (DipoleVector(0.0, 0.0, 1.0),)]
            @test_throws ArgumentError VectorPSF(1.33, 0.68, args...; n_medium=1.33,
                n_coverslip=1.33, n_immersion=1.33, grid_size=33)
            p = VectorPSF(1.33, 0.68, args...; n_medium=1.33,
                n_coverslip=1.52, n_immersion=1.52, grid_size=33)
            @test all(v -> all(isfinite, v.Ex.field) && all(isfinite, v.Ey.field),
                      p.vector_pupils)
            @test sum(energy, p.vector_pupils) ≈ 1
            @test all(x -> isfinite(p(x...)), points)
            bottom = isempty(args) ? p : VectorPSF(1.33, 0.68,
                DipoleVector(p.dipole.px, p.dipole.py, -p.dipole.pz);
                n_medium=1.33, n_coverslip=1.52, n_immersion=1.52, grid_size=33)
            @test isfinite(FourPiPSF(p, bottom)(points[2]...))
        end
    end

    @testset "I7 amplitude integration entry points" begin
        camera = IdealCamera(4, 3, 0.1)
        e = Emitter3D(0.2, 0.15, 0.1, 3.0)
        v = vector(DipoleVector(1.0, 0.0, 0.0))
        s = ScalarPSF(1.2, 0.68, 1.33)
        for p in (PlanePSF(v; dz=0.1), AnalyzerPSF(v; component=1),
                  FourPiPSF(v, v; phase=0.3), PlanePSF(s; dz=0.1),
                  FourPiPSF(s, s; phase=0.3))
            edges = (camera.pixel_edges_x, camera.pixel_edges_y)
            for threaded in (false, true)
                kw = (; sampling=2, threaded)
                expected = integrate_pixels_amplitude(p, camera, e; kw...)
                @test integrate_pixels_amplitude(p, edges..., e; kw...) == expected
                for geometry in ((camera,), edges)
                    out = similar(expected)
                    @test M.integrate_pixels_amplitude!(out, p, geometry..., e;
                        kw...) === out
                    @test out == expected
                    @test integrate_pixels_amplitude(p, geometry..., [e, e]; kw...) ≈
                          2 .* expected
                    @test size(integrate_pixels_amplitude(p, geometry...,
                        Emitter3D[]; kw...)) == size(expected)
                end
            end
        end
    end

    @testset "I8 default spline dimension trait" begin
        p = GaussianPSF(0.15)
        s = SplinePSF(p)
        @test s.z_range === nothing
        @test s(0.0, 0.0) ≈ p(0.0, 0.0)
        @test SplinePSF(s).z_range === nothing
    end

    @testset "I9 promoted vector precision" begin
        for args in [(), (DipoleVector(1f0, 0f0, 0f0),)]
            p = VectorPSF(1.2f0, 0.68f0, args...; n_medium=1.33f0,
                n_coverslip=1.33f0, n_immersion=1.33f0, grid_size=16)
            @test p isa VectorPSF{Float32}
            @test all(v -> eltype(v.Ex.field) == ComplexF32, p.vector_pupils)
            @test isfinite(p(0f0, 0f0, 0.1f0))
            c32 = ZernikeCoefficients(6; T=Float32)
            c32.phase[6] = 0.2f0
            kw = (; n_medium=1.33f0, n_coverslip=1.33f0,
                n_immersion=1.33f0, grid_size=16)
            aberrated = VectorPSF(1.2f0, 0.68f0, args...;
                kw..., base_zernike=c32)
            @test aberrated isa VectorPSF{Float32}
            @test isfinite(aberrated(0f0, 0f0, 0.1f0))
            mixed = VectorPSF(1.2, 0.68, args...; kw..., base_zernike=c32)
            @test mixed isa VectorPSF{Float64}
            @test mixed.zernike_coeffs isa ZernikeCoefficients{Float64}
            base32 = PupilFunction(1.2f0, 0.68f0, 1.33f0, c32; grid_size=16)
            mixed = VectorPSF(1.2, 0.68, args...; kw..., base_pupil=base32)
            @test mixed.base_pupil isa PupilFunction{Float64}
        end
        p = VectorPSF(1.2f0, 0.68f0, DipoleVector(1f0, 0f0, 0f0);
            n_medium=1.33f0, n_coverslip=1.33, n_immersion=1.33f0,
            z_stage=0f0, grid_size=16)
        @test p isa VectorPSF{Float64}
    end

    @testset "I10 coefficient assignment" begin
        zc = ZernikeCoefficients(15)
        zc[6, :phase] = 0.3
        zc[4, :mag] = 0.2
        @test zc.phase[6] == 0.3
        @test zc.mag[4] == 0.2
        @test_throws ArgumentError setindex!(zc, 0.1, 1, :unknown)
    end

    @testset "I10 README computational examples" begin
        readme = read(joinpath(@__DIR__, "..", "README.md"), String)
        blocks = collect(eachmatch(r"```julia\r?\n(.*?)```"s, readme))
        @test length(blocks) == 4
        for (i, block) in enumerate(blocks)
            code = block.captures[1]
            # Installation is exercised separately in an isolated package environment.
            occursin("Pkg.add", code) && continue
            scope = Module(gensym(:ReadmeExample))
            @test Base.include_string(scope, code, "README block $i") !==
                  nothing
        end
    end

    @testset "I11 finite support pixel bounds" begin
        edges = collect(0.0:0.1:1.0)
        e = Emitter2D(0.45, 0.45, 1.0)
        indices(s) = M.get_pixel_indices(edges, edges, e, s)
        @test indices((0.4, 0.5, 0.4, 0.5)) == (5:5, 5:5)
        @test indices((0.42, 0.48, 0.42, 0.48)) == (5:5, 5:5)
        @test indices(0.1) == (4:6, 4:6)
        @test indices(Inf) == (1:10, 1:10)
        @test all(isempty, indices((1.0, 1.2, 0.1, 0.2)))
        @test all(isempty, indices((-0.2, 0.0, 0.1, 0.2)))
        p = GaussianPSF(0.15)
        full = integrate_pixels(p, edges, edges, e)
        finite = integrate_pixels(p, edges, edges, e; support=(0.4, 0.5, 0.4, 0.5))
        @test count(!iszero, finite) == 1
        @test finite[5, 5] == full[5, 5]
    end

    @testset "I14 Gaussian width invariant" begin
        for bad in (0.0, -0.15, Inf, -Inf, NaN)
            @test_throws ArgumentError GaussianPSF(bad)
            @test_throws ArgumentError GaussianPSF{Float64}(bad)
        end
        @test GaussianPSF(0.15f0) isa GaussianPSF{Float32}
        @test GaussianPSF(0.15)(0.0, 0.0) == 1 / (2π * 0.15^2)
    end

    @testset "I15 vector tests cannot silently skip" begin
        for name in ("multi_emitter", "thread", "io")
            source = read(joinpath(@__DIR__, name * "_tests.jl"), String)
            @test !occursin(r"Skipping VectorPSF", source)
        end
    end

    @testset "I17 dependency declarations" begin
        project = read(joinpath(@__DIR__, "..", "Project.toml"), String)
        deps = split(split(project, "[deps]")[2], "[compat]")[1]
        compat = split(split(project, "[compat]")[2], "[extras]")[1]
        extras = split(split(project, "[extras]")[2], "[targets]")[1]
        @test !occursin("JLD2", project)
        @test !occursin("Statistics", deps)
        @test occursin("Statistics", extras)
        @test occursin("LinearAlgebra = \"1\"", compat)
        @test occursin("Test = \"1\"", compat)
        @test occursin(r"test = .*\"Statistics\"", project)
    end
end
