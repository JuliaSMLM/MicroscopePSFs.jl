using Test, MicroscopePSFs, HDF5
import Base.CoreLogging as Logging

part4_fields(p::VectorPSF) = [copy(f.field)
    for v in p.vector_pupils for f in (v.Ex, v.Ey)]
part4_fields(p::Union{PlanePSF,AnalyzerPSF}) = part4_fields(p.psf)
part4_fields(p::FourPiPSF) = vcat(part4_fields(p.top), part4_fields(p.bottom))
part4_fields(p::MultiChannelPSF) = reduce(vcat, part4_fields.(p.channels))
part4_fields(p::VectorPupilFunction) = [copy(p.Ex.field), copy(p.Ey.field)]

@testset "0.5.7 pupil operations and typed dipoles" verbose=true begin
    M = MicroscopePSFs
    fixture(name) = joinpath(@__DIR__, "fixtures", "vector-v056-$name.h5")
    vector(args...) = VectorPSF(1.2, 0.68, args...; grid_size=16,
        n_medium=1.33, n_coverslip=1.33, n_immersion=1.33,
        base_zernike=ZernikeCoefficients(6))

    @testset "0.5.6 fixtures load silently and preserve stored pupil fields" begin
        for name in ("tilted", "bottom", "rotating", "transverse", "standalone")
            p = @test_logs min_level=Logging.Warn load_psf(fixture(name))
            h5open(fixture(name), "r") do f
                if p isa VectorPupilFunction
                    @test p.Ex.field == M._io_load_complex_array(f["ex_field"], "field")
                    @test p.Ey.field == M._io_load_complex_array(f["ey_field"], "field")
                else
                    for (i, pupil) in enumerate(p.vector_pupils)
                        group = f["data/pupil_$i"]
                        @test pupil.Ex.field ==
                            M._io_load_complex_array(group["pupil_ex"], "field")
                        @test pupil.Ey.field ==
                            M._io_load_complex_array(group["pupil_ey"], "field")
                    end
                end
            end
        end
    end

    @testset "0.5.6 fields survive repeated saves and nested containers" begin
        # Fixtures were written in a separate process from git archive v0.5.6.
        # See fixtures/write_v056.jl for the reproducible writer command.
        old = @test_logs min_level=Logging.Warn load_psf(fixture("tilted"))
        bottom = @test_logs min_level=Logging.Warn load_psf(fixture("bottom"))
        rotating = @test_logs min_level=Logging.Warn load_psf(fixture("rotating"))
        transverse = @test_logs min_level=Logging.Warn load_psf(fixture("transverse"))
        current = vector(DipoleVector(0.6, 0., 0.8))
        current_bottom = vector(DipoleVector(0.6, 0., -0.8))
        matched = @test_logs min_level=Logging.Warn FourPiPSF(old, bottom)
        mixed = @test_logs min_level=Logging.Warn FourPiPSF(
            PlanePSF(old; dz=0.1), AnalyzerPSF(current_bottom; component=1))
        containers = (old, PlanePSF(old; dz=0.1), AnalyzerPSF(old; component=1),
            matched, mixed, rotating,
            MultiChannelPSF(; channels=[PlanePSF(old; dz=0.1), current]),
            MultiChannelPSF(; channels=[matched, PlanePSF(rotating; dz=0.2), current]))
        mktempdir() do dir
            for (i, p) in enumerate(containers)
                path = joinpath(dir, "nested-$i.h5")
                original_fields = part4_fields(p)
                save_psf(path, p)
                loaded = @test_logs min_level=Logging.Warn load_psf(path)
                @test part4_fields(loaded) == original_fields
                h5open(path, "r") do f
                    @test read(attributes(f)["io_version"]) == "1.1.0"
                end
                if !(p isa MultiChannelPSF)
                    @test loaded(0.1, -0.2, 0.3) == p(0.1, -0.2, 0.3)
                end
                # A second save must preserve the exact fields.
                save_psf(path, loaded)
                again = @test_logs min_level=Logging.Warn load_psf(path)
                @test part4_fields(again) == original_fields
                # Check nested field loading with the earlier supported file version.
                if p isa Union{PlanePSF,FourPiPSF,MultiChannelPSF}
                    h5open(path, "r+") do f
                        HDF5.delete_attribute(f, "io_version")
                        attributes(f)["io_version"] = "1.0.0"
                    end
                    unmarked = @test_logs min_level=Logging.Warn load_psf(path)
                    @test part4_fields(unmarked) == original_fields
                end
            end
            for p in (current, transverse, PlanePSF(current; dz=0.1),
                FourPiPSF(current, current_bottom),
                MultiChannelPSF(; channels=[current, transverse]))
                path = joinpath(dir, "current.h5")
                save_psf(path, p)
                loaded = @test_logs min_level=Logging.Warn load_psf(path)
                @test part4_fields(loaded) == part4_fields(p)
            end
            path = joinpath(dir, "unversioned.h5")
            cp(fixture("tilted"), path)
            h5open(path, "r+") do f
                HDF5.delete_attribute(f, "io_version")
            end
            loaded = @test_logs min_level=Logging.Warn load_psf(path)
            @test part4_fields(loaded) == part4_fields(old)
        end
        @test part4_fields(old) != part4_fields(current)
        @test h5open(f -> read(attributes(f)["io_version"]), fixture("tilted")) ==
            "1.0.0"
    end

    @testset "Pupil construction, copying, defocus and rebuilding" begin
        for name in ("tilted", "rotating")
            old = @test_logs min_level=Logging.Warn load_psf(fixture(name))
            original = part4_fields(old)
            for vp in old.vector_pupils
                copied = copy(vp)
                @test part4_fields(copied) == part4_fields(vp)
                @test copied.Ex.field !== vp.Ex.field
                @test copied.Ey.field !== vp.Ey.field
                @test M.apply_defocus!(copied, 0.2, 0.1) === copied
                default = VectorPupilFunction(vp.nₐ, vp.λ, vp.n_medium,
                    vp.n_coverslip, vp.n_immersion, copy(vp.Ex), copy(vp.Ey))
                @test part4_fields(default) == part4_fields(vp)
                mktempdir() do dir
                    path = joinpath(dir, "pupil.h5")
                    save_psf(path, default)
                    loaded = @test_logs min_level=Logging.Warn load_psf(path)
                    @test part4_fields(loaded) == part4_fields(default)
                    M.fill_vector_pupils!(loaded, DipoleVector(0.6, 0., 0.8))
                    save_psf(path, loaded)
                    @test_logs min_level=Logging.Warn load_psf(path)
                end
            end
            pupil_at(old, 0.2; z_stage=0.1)
            @test part4_fields(old) == original
            @test M.update_pupils!(old) === old
            rebuilt = name == "rotating" ? vector() : vector(old.dipole)
            @test part4_fields(old) ≈ part4_fields(rebuilt) rtol=1e-14
            mktempdir() do dir
                path = joinpath(dir, "rebuilt.h5")
                save_psf(path, old)
                @test_logs min_level=Logging.Warn load_psf(path)
            end
        end
    end

    @testset "Explicit dipole and emitter types are retained" begin
        components(d) = (d.px, d.py, d.pz)
        for T in (Int8, UInt8, Int, BigInt), args in ((1,0,0), (0,2,0), (0,0,3), (0,0,0))
            d = DipoleVector{T}(args...)
            @test d isa DipoleVector{T}
            @test eltype(DipoleVector{T}[d]) == DipoleVector{T}
            @test components(d) == map(a -> iszero(a) ? 0 : 1, args)
            if any(!iszero, args)
                e = DipoleEmitter3D{T}(1, 2, 0, 3, d)
                @test e isa DipoleEmitter3D{T}
                @test e.dipole isa DipoleVector{T}
                @test (e.x, e.y, e.z, e.photons) == (1, 2, 0, 3)
                @test eltype(DipoleEmitter3D{T}[e]) == DipoleEmitter3D{T}
            else
                @test_throws ArgumentError DipoleEmitter3D{T}(1, 2, 0, 3, d)
            end
        end
        @test components(DipoleVector{Int}(-2,0,0)) == (-1,0,0)
        @test components(DipoleVector{Int}(0.,0.25,0.)) == (0,1,0)
        @test components(DipoleVector{Int8}(300,0,0)) == (1,0,0)
        @test components(DipoleVector{Int}(floatmax(Float64),0,0)) == (1,0,0)
        for T in (Int8, UInt8, Int, BigInt), args in ((1,1,0), (2,3,6), (3,4,0))
            @test_throws r"floating-point type.*untyped DipoleVector" DipoleVector{T}(
                args...)
            @test_throws r"floating-point type.*untyped DipoleEmitter3D" DipoleEmitter3D{T}(
                1, 2, 0, 3, DipoleVector(args...))
        end
        for T in (Rational{Int}, Rational{BigInt})
            @test components(DipoleVector{T}(3,4,0)) == (3//5,4//5,0//1)
            @test DipoleVector{T}(0,0,0) isa DipoleVector{T}
            @test_throws ArgumentError DipoleVector{T}(1,1,0)
            e = DipoleEmitter3D{T}(1//2, 2, 0, 3, DipoleVector{T}(3,4,0))
            @test e isa DipoleEmitter3D{T}
            @test components(e.dipole) == (3//5,4//5,0//1)
        end
        @test_throws ArgumentError DipoleVector{UInt8}(-1,0,0)
        @test_throws ArgumentError DipoleEmitter3D{Int}(0.5, 0, 0, 1,
            DipoleVector(1,0,0))
        for args in ((2,3,6), (2//1,3//1,6//1), (0,0,0))
            @test DipoleVector(args...) isa DipoleVector{Float64}
        end
        for T in (Float32, Float64, BigFloat)
            d = DipoleVector{T}(T(2),T(3),T(6))
            @test d isa DipoleVector{T}
            @test sum(abs2, components(d)) ≈ 1
            @test DipoleEmitter3D{T}(1, 2, 0, 3, d) isa DipoleEmitter3D{T}
        end
        @test DipoleEmitter3D(1, 2, 0, 3, DipoleVector{Int}(1,0,0)) isa
            DipoleEmitter3D{Float64}
    end
end
