using Test, MicroscopePSFs, HDF5

struct CachedReviewPSF <: MicroscopePSFs.Abstract3DPSF{Float64}
    fields::Vector{Vector{ComplexF64}}
end
MicroscopePSFs._n_components(::CachedReviewPSF) = 2
MicroscopePSFs.dipole_fields(p::CachedReviewPSF, x::Real, y::Real, z::Real) = p.fields

@testset "0.5.7 rulings and ecosystem regressions" verbose=true begin
    M = MicroscopePSFs
    warn = Base.CoreLogging.Warn
    vector(args...) = VectorPSF(1.2, 0.68, args...; grid_size=16,
        n_medium=1.33, n_coverslip=1.33, n_immersion=1.33)
    scalar = ScalarPSF(1.2, 0.68, 1.33)

    @testset "I6 file version identifies the axial sign" begin
        @test M.IO_VERSION == v"1.1.0"
        @test M.MIN_SUPPORTED_VERSION == v"1.0.0"
        transverse = vector(DipoleVector(1.0, 0.0, 0.0))
        axial = vector(DipoleVector(0.0, 0.0, 1.0))
        rotating = vector()
        cavity = FourPiPSF(rotating, rotating)
        channels = MultiChannelPSF(; channels=[PlanePSF(axial; dz=0.1), cavity])
        fixtures = ((GaussianPSF(0.15), false), (scalar, false),
            (transverse, false), (axial, true), (rotating, true),
            (PlanePSF(AnalyzerPSF(cavity; component=1); dz=0.1), true),
            (channels, true), (MultiChannelPSF(; channels=[scalar, transverse]), false))
        function check_groups(group)
            @test !haskey(attributes(group), "axial_sign_convention")
            for name in keys(group)
                child = group[name]
                child isa HDF5.Group && check_groups(child)
            end
        end
        fields(::Any) = Matrix{ComplexF64}[]
        fields(p::VectorPSF) = [copy(f.field)
            for v in p.vector_pupils for f in (v.Ex, v.Ey)]
        fields(p::Union{PlanePSF,AnalyzerPSF}) = fields(p.psf)
        fields(p::FourPiPSF) = vcat(fields(p.top), fields(p.bottom))
        fields(p::MultiChannelPSF) = reduce(vcat, fields.(p.channels))
        message = Regex(raw"no axial-sign marker.*0\.5\.6 or earlier.*" *
            raw"old axial dipole sign.*unreleased 0\.5\.7 development build.*" *
            raw"cannot be told from the file.*fields are preserved.*" *
            raw"rebuild.*parameters to be sure")
        for (p, affected) in fixtures
            mktempdir() do dir
                path = joinpath(dir, "convention.h5")
                save_psf(path, p)
                h5open(path, "r") do f
                    @test read(attributes(f)["io_version"]) == "1.1.0"
                    check_groups(f)
                end
                loaded = @test_logs min_level=warn load_psf(path)
                @test fields(loaded) == fields(p)
                h5open(path, "r+") do f
                    HDF5.delete_attribute(f, "io_version")
                    attributes(f)["io_version"] = "1.0.0"
                end
                loaded = if affected
                    @test_logs (:warn, message) min_level=warn load_psf(path)
                else
                    @test_logs min_level=warn load_psf(path)
                end
                @test fields(loaded) == fields(p)
                if affected
                    h5open(path, "r+") do f
                        HDF5.delete_attribute(f, "io_version")
                    end
                    @test_logs (:warn, message) min_level=warn load_psf(path)
                end
            end
        end
        @test_logs min_level=warn M._check_version_compatibility("1.1.0")
    end

    @testset "A1 every dipole constructor normalizes" begin
        components(d) = [d.px, d.py, d.pz]
        for args in ((2.0, 0.0, 0.0), (2f0, 3f0, 6f0), (1.0, 0, 1.0),
                     (2, 3, 6), (2//1, 3//1, 6//1))
            d = DipoleVector(args...)
            @test sum(abs2, components(d)) ≈ 1
            @test components(d) ≈ collect(args) ./ sqrt(sum(abs2, args))
        end
        for T in (Float32, Float64)
            d = DipoleVector{T}(2, 3, 6)
            @test d isa DipoleVector{T}
            @test components(d) ≈ T[2, 3, 6] ./ 7
            z = DipoleVector{T}(0, 0, 0)
            @test all(iszero, components(z))
            @test_throws ArgumentError vector(z)
            @test_throws ArgumentError DipoleEmitter3D(1.0, 1.0, 0.0, 1.0, z)
            @test_throws ArgumentError DipoleEmitter3D{T}(1, 1, 0, 1, z)
        end
        @test_throws ArgumentError DipoleEmitter3D(1., 1., 0., 1., 0., 0., 0.)
        @test DipoleEmitter3D(1., 1., 0., 1., 2., 0., 0.).dipole.px == 1
        @test DipoleEmitter3D(1., 1., 0., 1., DipoleVector(2., 0., 0.)).dipole.px == 1
        @test all(iszero, components(DipoleVector(0, 0f0, 0.0)))
        @test all(iszero, components(vector().dipole))
        @test vector(DipoleVector(2., 0., 0.))(0.1, 0.0, 0.2) ≈
              vector(DipoleVector(1., 0., 0.))(0.1, 0.0, 0.2)
        for scale in (floatmax(Float64), floatmin(Float64))
            @test sum(abs2, components(DipoleVector(scale, scale, scale))) ≈ 1
        end
    end

    @testset "A2 unit-photon amplitude documentation" begin
        doc = lowercase(string(@doc integrate_pixels_amplitude))
        @test occursin("unit-photon", doc)
        @test occursin("photons are ignored", doc)
        camera = IdealCamera(2, 2, 0.1)
        e1 = Emitter2D(0.1, 0.1, 1.0)
        e4 = Emitter2D(0.1, 0.1, 4.0)
        psf = GaussianPSF(0.15)
        a = integrate_pixels_amplitude(psf, camera, e1)
        @test a == integrate_pixels_amplitude(psf, camera, e4)
        @test !isapprox(abs2.(a), integrate_pixels(psf, camera, e1))
    end

    @testset "A3 unused axial-phase parameters documented" begin
        doc = lowercase(string(@doc calculate_axial_phase))
        @test occursin("ignored", doc)
        @test occursin("kz_coverslip", doc) && occursin("coverslip_thickness", doc)
        a = calculate_axial_phase(0.3, 0.1, 1.9+0im, 2.0+0im, 2.2+0im, 0.17)
        b = calculate_axial_phase(0.3, 0.1, 1.9+0im, 9.0+0im, 2.2+0im, 5.0)
        @test a == b
    end

    @testset "4Pi photon and axial-offset documentation" begin
        for doc in (string(@doc MultiChannelPSF), string(@doc channel_position),
                    read(joinpath(@__DIR__, "..", "docs", "src", "psfs",
                        "multichannel.md"), String))
            text = lowercase(replace(doc, r"\s+" => " "))
            @test occursin("fringe-average", text)
            @test occursin("twice", text)
            @test occursin("before those losses", text)
            @test occursin("learned", text) && occursin("z_offset", text)
        end
    end

    @testset "Numeric channel collections accept tuples" begin
        reg = ([1f0 0 0; 0 1 0], [1f0 0 0.2f0; 0 1 0.3f0])
        mc = MultiChannelPSF(; channels=(scalar, scalar), registration=reg,
            transmission=(0.5f0, 0.5f0), z_offset=(0f0, 0.2f0))
        @test mc isa MultiChannelPSF{Float32}
        @test mc.transmission == Float32[0.5, 0.5]
        @test mc.z_offset == Float32[0, 0.2]
        @test mc.registration == collect(reg)
        @test channel_position(mc, 2, 0f0, 0f0, 0f0) == (0.2f0, 0.3f0, 0.2f0)
    end

    @testset "Transmission tolerance covers summation rounding" begin
        for T in (Float32, Float64)
            @test_throws ArgumentError MultiChannelPSF(; channels=[scalar],
                transmission=T[1 + sqrt(eps(T))/2])
            for N in (1, 2, 10, 100)
                channels = fill(scalar, N)
                mc = MultiChannelPSF(; channels, transmission=fill(one(T)/N, N))
                @test mc isa MultiChannelPSF{T}
                almost = vcat(T[1 + 4N*eps(T)], zeros(T, N-1))
                @test MultiChannelPSF(; channels, transmission=almost) isa MultiChannelPSF
                excess = vcat(T[1 + 16N*eps(T)], zeros(T, N-1))
                @test_throws ArgumentError MultiChannelPSF(; channels, transmission=excess)
            end
        end
    end

    @testset "Analyzer fields leave cached inner arrays intact" begin
        original = [ComplexF64[1+2im, 3+4im], ComplexF64[5+6im, 7+8im]]
        cached = CachedReviewPSF(deepcopy(original))
        h, v = AnalyzerPSF(cached; component=1), AnalyzerPSF(cached; component=2)
        horizontal = dipole_fields(h, 0., 0., 0.)
        @test cached.fields == original
        vertical = dipole_fields(v, 0., 0., 0.)
        @test cached.fields == original
        @test horizontal == [[e[1], 0] for e in original]
        @test vertical == [[0, e[2]] for e in original]
        horizontal[1][1] = 0
        @test cached.fields == original
        @test vertical[1][2] == original[1][2]
    end
end
