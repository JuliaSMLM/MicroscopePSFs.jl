using Test, MicroscopePSFs, HDF5

part4_fields(p::VectorPSF) = [copy(f.field)
    for v in p.vector_pupils for f in (v.Ex, v.Ey)]
part4_fields(p::Union{PlanePSF,AnalyzerPSF}) = part4_fields(p.psf)
part4_fields(p::FourPiPSF) = vcat(part4_fields(p.top), part4_fields(p.bottom))
part4_fields(p::MultiChannelPSF) = reduce(vcat, part4_fields.(p.channels))
part4_fields(p::VectorPupilFunction) = [copy(p.Ex.field), copy(p.Ey.field)]
part4_conventions(p::VectorPSF) = part4_conventions.(p.vector_pupils)
part4_conventions(p::Union{PlanePSF,AnalyzerPSF}) = part4_conventions(p.psf)
part4_conventions(p::FourPiPSF) =
    vcat(part4_conventions(p.top), part4_conventions(p.bottom))
part4_conventions(p::MultiChannelPSF) = reduce(vcat, part4_conventions.(p.channels))
part4_conventions(p::VectorPupilFunction) =
    hasproperty(p, :axial_sign_convention) ? p.axial_sign_convention : :missing

function part4_stored_conventions(group)
    attrs = attributes(group)
    result = haskey(attrs, "axial_sign_convention") ?
        [read(attrs["axial_sign_convention"])] : String[]
    for name in keys(group)
        child = group[name]
        child isa HDF5.Group && append!(result, part4_stored_conventions(child))
    end
    return result
end

function part4_erase_markers!(group)
    attrs = attributes(group)
    if haskey(attrs, "axial_sign_convention")
        HDF5.delete_attribute(group, "axial_sign_convention")
    end
    for name in keys(group)
        child = group[name]
        child isa HDF5.Group && part4_erase_markers!(child)
    end
    return nothing
end

@testset "0.5.7 pupil provenance and typed dipoles" verbose=true begin
    M = MicroscopePSFs
    warn = Base.CoreLogging.Warn
    message = Regex(raw"no axial-sign marker.*0\.5\.6 or earlier.*" *
        raw"old axial dipole sign.*unreleased 0\.5\.7 development build.*" *
        raw"cannot be told from the file.*fields are preserved.*" *
        raw"rebuild.*parameters to be sure")
    mismatch = r"4Pi arms have different axial sign conventions.*interference.*Rebuild"
    fixture(name) = joinpath(@__DIR__, "fixtures", "vector-v056-$name.h5")
    vector(args...) = VectorPSF(1.2, 0.68, args...; grid_size=16,
        n_medium=1.33, n_coverslip=1.33, n_immersion=1.33,
        base_zernike=ZernikeCoefficients(6))

    @testset "Actual 0.5.6 fields stay marked after re-saving and nesting" begin
        # Fixtures were written in a separate process from git archive v0.5.6.
        # See fixtures/write_v056.jl for the reproducible writer command.
        old = @test_logs (:warn, message) min_level=warn load_psf(fixture("tilted"))
        bottom = @test_logs (:warn, message) min_level=warn load_psf(fixture("bottom"))
        rotating = @test_logs (:warn, message) min_level=warn load_psf(fixture("rotating"))
        transverse = @test_logs min_level=warn load_psf(fixture("transverse"))
        current = vector(DipoleVector(0.6, 0., 0.8))
        current_bottom = vector(DipoleVector(0.6, 0., -0.8))
        matched = @test_logs min_level=warn FourPiPSF(old, bottom)
        mixed = @test_logs (:warn, mismatch) min_level=warn FourPiPSF(
            PlanePSF(old; dz=0.1), AnalyzerPSF(current_bottom; component=1))
        containers = (old, PlanePSF(old; dz=0.1), AnalyzerPSF(old; component=1),
            matched, mixed, rotating,
            MultiChannelPSF(; channels=[PlanePSF(old; dz=0.1), current]),
            MultiChannelPSF(; channels=[matched, PlanePSF(rotating; dz=0.2), current]))
        mktempdir() do dir
            for (i, p) in enumerate(containers)
                path = joinpath(dir, "legacy-$i.h5")
                original_fields = part4_fields(p)
                save_psf(path, p)
                loaded = @test_logs (:warn, message) min_level=warn load_psf(path)
                @test part4_fields(loaded) == original_fields
                @test part4_conventions(loaded) == part4_conventions(p)
                h5open(path, "r") do f
                    @test read(attributes(f)["io_version"]) == "1.1.0"
                    @test sort(part4_stored_conventions(f)) ==
                        sort(string.(part4_conventions(p)))
                end
                if !(p isa MultiChannelPSF)
                    @test loaded(0.1, -0.2, 0.3) == p(0.1, -0.2, 0.3)
                end
                # A second re-save must preserve both the warning and exact fields.
                save_psf(path, loaded)
                again = @test_logs (:warn, message) min_level=warn load_psf(path)
                @test part4_fields(again) == original_fields
                # Exercise unmarked nested groups using these actual legacy fields.
                if p isa Union{PlanePSF,FourPiPSF,MultiChannelPSF}
                    h5open(path, "r+") do f
                        part4_erase_markers!(f)
                        HDF5.delete_attribute(f, "io_version")
                        attributes(f)["io_version"] = "1.0.0"
                    end
                    unmarked = @test_logs (:warn, message) min_level=warn load_psf(path)
                    @test part4_fields(unmarked) == original_fields
                    @test all(==(:legacy), part4_conventions(unmarked))
                end
            end
            for p in (current, transverse, PlanePSF(current; dz=0.1),
                FourPiPSF(current, current_bottom),
                MultiChannelPSF(; channels=[current, transverse]))
                path = joinpath(dir, "current.h5")
                save_psf(path, p)
                loaded = @test_logs min_level=warn load_psf(path)
                @test part4_fields(loaded) == part4_fields(p)
                @test all(==(:current), part4_conventions(loaded))
            end
            # Missing versions follow the same conservative rule as format 1.0.0.
            path = joinpath(dir, "unversioned.h5")
            cp(fixture("tilted"), path)
            h5open(path, "r+") do f
                HDF5.delete_attribute(f, "io_version")
            end
            loaded = @test_logs (:warn, message) min_level=warn load_psf(path)
            @test part4_conventions(loaded) == [:legacy]
        end
        @test part4_conventions(old) == [:legacy]
        @test part4_conventions(rotating) == fill(:legacy, 3)
        @test part4_conventions(transverse) == [:current]
        @test part4_conventions(current) == [:current]
        @test part4_fields(old) != part4_fields(current)
        @test h5open(f -> read(attributes(f)["io_version"]), fixture("tilted")) ==
            "1.0.0"
    end

    @testset "Pupil construction, copying, defocus and rebuilding" begin
        for name in ("tilted", "rotating")
            old = @test_logs (:warn, message) min_level=warn load_psf(fixture(name))
            original = part4_fields(old)
            for vp in old.vector_pupils
                copied = copy(vp)
                @test part4_conventions(copied) == :legacy
                @test part4_fields(copied) == part4_fields(vp)
                @test copied.Ex.field !== vp.Ex.field
                @test copied.Ey.field !== vp.Ey.field
                @test M.apply_defocus!(copied, 0.2, 0.1) === copied
                @test part4_conventions(copied) == :legacy
                default = VectorPupilFunction(vp.nₐ, vp.λ, vp.n_medium,
                    vp.n_coverslip, vp.n_immersion, copy(vp.Ex), copy(vp.Ey))
                @test part4_conventions(default) == :current
                marked = VectorPupilFunction(vp.nₐ, vp.λ, vp.n_medium,
                    vp.n_coverslip, vp.n_immersion, copy(vp.Ex), copy(vp.Ey);
                    axial_sign_convention=:legacy)
                @test part4_fields(marked) == part4_fields(vp)
                @test_throws ArgumentError VectorPupilFunction(vp.nₐ, vp.λ,
                    vp.n_medium, vp.n_coverslip, vp.n_immersion, vp.Ex, vp.Ey;
                    axial_sign_convention=:unknown)
                mktempdir() do dir
                    path = joinpath(dir, "pupil.h5")
                    save_psf(path, marked)
                    loaded = @test_logs (:warn, message) min_level=warn load_psf(path)
                    @test part4_conventions(loaded) == :legacy
                    @test part4_fields(loaded) == part4_fields(marked)
                    M.fill_vector_pupils!(loaded, DipoleVector(0.6, 0., 0.8))
                    @test part4_conventions(loaded) == :current
                    save_psf(path, loaded)
                    @test_logs min_level=warn load_psf(path)
                end
            end
            defocused = pupil_at(old, 0.2; z_stage=0.1)
            pupils = defocused isa VectorPupilFunction ? [defocused] : defocused
            @test all(==(:legacy), part4_conventions.(pupils))
            @test part4_fields(old) == original
            @test M.update_pupils!(old) === old
            @test all(==(:current), part4_conventions(old))
            rebuilt = name == "rotating" ? vector() : vector(old.dipole)
            @test part4_fields(old) ≈ part4_fields(rebuilt) rtol=1e-14
            mktempdir() do dir
                path = joinpath(dir, "rebuilt.h5")
                save_psf(path, old)
                @test_logs min_level=warn load_psf(path)
            end
        end
        a, b = vector(), vector()
        b.vector_pupils[3].axial_sign_convention = :legacy
        @test_logs (:warn, mismatch) min_level=warn FourPiPSF(a, b)
        @test_logs (:warn, mismatch) min_level=warn FourPiPSF(
            AnalyzerPSF(a; component=1), PlanePSF(b; dz=0.1))
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
