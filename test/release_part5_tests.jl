using Test, MicroscopePSFs, HDF5
import Base.CoreLogging as Logging
include("angular_reference.jl")

part5_fields(p::VectorPupilFunction) = (copy(p.Ex.field), copy(p.Ey.field))

@testset "PR 40 second-review regressions" verbose=true begin
    M = MicroscopePSFs
    fixture(name) = joinpath(@__DIR__, "fixtures", "vector-v056-$name.h5")

    @testset "Standalone 0.5.6 pupil fields survive repeated saves" begin
        p = @test_logs min_level=Logging.Warn load_psf(fixture("standalone"))
        original = part5_fields(p)
        h5open(fixture("standalone"), "r") do f
            @test read(attributes(f)["io_version"]) == "1.0.0"
            @test p.Ex.field == M._io_load_complex_array(f["ex_field"], "field")
            @test p.Ey.field == M._io_load_complex_array(f["ey_field"], "field")
        end
        mktempdir() do dir
            path = joinpath(dir, "resaved.h5")
            for _ in 1:2
                save_psf(path, p)
                h5open(path, "r") do f
                    @test read(attributes(f)["io_version"]) == "1.1.0"
                end
                p = @test_logs min_level=Logging.Warn load_psf(path)
                @test part5_fields(p) == original
            end
            cp(fixture("standalone"), path; force=true)
            h5open(path, "r+") do f
                HDF5.delete_attribute(f, "io_version")
            end
            q = @test_logs min_level=Logging.Warn load_psf(path)
            @test part5_fields(q) == original
            M.fill_vector_pupils!(q, DipoleVector(0.6,0.,0.8))
            save_psf(path, q)
            current = @test_logs min_level=Logging.Warn load_psf(path)
            @test part5_fields(current) == part5_fields(q)
            h5open(path, "r+") do f
                HDF5.delete_attribute(f, "io_version")
                attributes(f)["io_version"] = "1.0.0"
            end
            reloaded = @test_logs min_level=Logging.Warn load_psf(path)
            @test part5_fields(reloaded) == part5_fields(current)
        end
    end

    @testset "Matched-index domain and angular-reference accuracy" begin
        points = [(0.,0.,0.), (0.1,-0.2,0.3), (0.2,0.1,0.), (0.,0.,0.3)]
        n, λ = 1.518, 0.68
        limit = 0.995*n
        for args in ((), (DipoleVector(0.6,0.,0.8),))
            for N in (32,33,64,65,128,129,255,256)
                for NA in (prevfloat(limit),limit)
                    p = VectorPSF(NA,λ,args...; n_medium=n,n_coverslip=n,
                        n_immersion=n,grid_size=N)
                    @test all(isfinite, [p(point...) for point in points])
                    energy = sum(M._vector_pupil_energy, p.vector_pupils)
                    @test energy ≈ 1 atol=4e-13
                end
                error = try
                    VectorPSF(nextfloat(limit),λ,args...; n_medium=n,n_coverslip=n,
                        n_immersion=n,grid_size=N)
                catch caught
                    caught
                end
                @test error isa ArgumentError
                if error isa ArgumentError
                    for text in ("0.995", "apodization diverges", "index-matched sample",
                        "grid sampling fails near 90 degrees", "real objectives")
                        @test occursin(text, error.msg)
                    end
                end
            end
        end
        vp = VectorPupilFunction(nextfloat(limit),λ,n,n,n,33)
        for normalize in (false,true)
            @test_throws ArgumentError M.fill_vector_pupils!(vp,DipoleVector(1,0,0);
                normalize=normalize)
        end
        # The rejection also catches the original catastrophic near-grazing case.
        @test_throws ArgumentError VectorPSF(prevfloat(1.33),λ; n_medium=1.33,
            n_coverslip=1.33,n_immersion=1.33,grid_size=33)
        for T in (Float32,Float64), N in (33,64)
            nm = T(n)
            boundary = T(0.995)*nm
            p = VectorPSF(boundary,T(λ); n_medium=nm,n_coverslip=nm,
                n_immersion=nm,grid_size=N)
            @test isfinite(p(0.1,-0.2,0.3))
            @test_throws ArgumentError VectorPSF(nextfloat(boundary),T(λ);
                n_medium=nm,n_coverslip=nm,n_immersion=nm,grid_size=N)
        end
        reference_base = part5_angular_reference(0.982*n,n,λ,points)
        reference_bound = part5_angular_reference(limit,n,λ,points)
        for (args,base_ref,bound_ref) in zip(((),(DipoleVector(0.6,0.,0.8),)),
            reference_base,reference_bound)
            for N in (64,65,128,129,255,256)
                a = VectorPSF(0.982*n,λ,args...; n_medium=n,n_coverslip=n,
                    n_immersion=n,grid_size=N)
                b = VectorPSF(limit,λ,args...; n_medium=n,n_coverslip=n,
                    n_immersion=n,grid_size=N)
                base_error = maximum(abs.([a(p...) for p in points]./base_ref .- 1))
                bound_error = maximum(abs.([b(p...) for p in points]./bound_ref .- 1))
                @test bound_error <= 2base_error
            end
        end
        # The real oil objective remains in the supported domain.
        oil = VectorPSF(1.49,λ; n_medium=n,n_coverslip=n,n_immersion=n)
        @test isfinite(oil(0.1,-0.2,0.3))
    end

    @testset "Cavity evaluation stays stable and observes refilled arms" begin
        top = @test_logs min_level=Logging.Warn load_psf(fixture("tilted"))
        bottom = @test_logs min_level=Logging.Warn load_psf(fixture("bottom"))
        cavity = @test_logs min_level=Logging.Warn FourPiPSF(top,bottom)
        M.update_pupils!(top)
        intensities = @test_logs min_level=Logging.Warn begin
            values = [cavity(0.1,-0.2,0.3) for _ in 1:3]
            amplitude(cavity,0.1,-0.2,0.3)
            dipole_fields(cavity,0.1,-0.2,0.3)
            PlanePSF(cavity; dz=0.1)(0.1,-0.2,0.2)
            AnalyzerPSF(cavity; component=1)(0.1,-0.2,0.3)
            values
        end
        @test all(==(first(intensities)), intensities)
        @test_logs min_level=Logging.Warn FourPiPSF(top,bottom)
        M.update_pupils!(bottom)
        @test_logs min_level=Logging.Warn cavity(0.1,-0.2,0.3)
        @test cavity(0.1,-0.2,0.3) != first(intensities)
        a = @test_logs min_level=Logging.Warn load_psf(fixture("rotating"))
        b = @test_logs min_level=Logging.Warn load_psf(fixture("rotating"))
        rotating = FourPiPSF(PlanePSF(a; dz=0.1), AnalyzerPSF(b; component=1))
        M.update_pupils!(b)
        @test_logs min_level=Logging.Warn rotating(0.1,-0.2,0.3)
    end

    @testset "Narrowing preserves finite nonzero dipole directions" begin
        components(d) = (d.px,d.py,d.pz)
        cases = ((Float32,(1e-50,2e-50,0.)), (Float16,(1e-8,3e-8,0.)),
            (Float32,(floatmax(Float64),0.,0.)),
            (Float16,(floatmax(Float64),floatmax(Float64),0.)),
            (Float32,(nextfloat(0.),nextfloat(0.),nextfloat(0.))))
        for (T,args) in cases
            d = DipoleVector{T}(args...)
            ref = DipoleVector(args...)
            @test d isa DipoleVector{T}
            @test components(d) === T.(components(ref))
            @test all(isfinite, components(d))
            @test any(!iszero, components(d))
            @test sum(abs2, Float64.(components(d))) ≈ 1 rtol=2eps(T)
            e = DipoleEmitter3D{T}(1,2,0,3,ref)
            @test e isa DipoleEmitter3D{T}
            @test components(e.dipole) === components(d)
        end
        @test components(DipoleVector{Float32}(1e-50,2e-50,0.)) ==
            (Float32(1/sqrt(5)),Float32(2/sqrt(5)),0f0)
        for T in (Float16,Float32,Float64), bad in (NaN,Inf,-Inf)
            @test_throws ArgumentError DipoleVector{T}(bad,0,0)
            @test all(iszero, components(DipoleVector{T}(0,0,0)))
        end
        for T in (Int8,Int), args in ((1,0,0),(0,0,0))
            @test DipoleVector{T}(args...) isa DipoleVector{T}
        end
        @test_throws ArgumentError DipoleVector{Int}(1,2,0)
        @test DipoleVector{Rational{Int}}(3,4,0).px == 3//5
        wide = DipoleVector{Float64}(2f0,3f0,6f0)
        @test components(wide) === components(DipoleVector(2.0,3.0,6.0))
    end
end
