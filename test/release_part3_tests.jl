using Test, MicroscopePSFs

function review_vector_energy(p)
    total = 0.0
    for pupil in p.vector_pupils
        N = size(pupil.Ex.field, 1)
        dk = 2 * p.nₐ / p.λ / (N - 1)
        for i in 1:N, j in 1:N
            kx, ky = (i - (N+1)/2) * dk, (j - (N+1)/2) * dk
            if (2i - N - 1)^2 + (2j - N - 1)^2 <= (N - 1)^2    # exact aperture geometry
                total += (abs2(pupil.Ex.field[j,i]) + abs2(pupil.Ey.field[j,i])) * dk^2
            end
        end
    end
    return total
end

@testset "PR 40 review regressions" verbose=true begin
    M = MicroscopePSFs
    fixed = (DipoleVector(1.0, 0.0, 1.0),)
    optics = (; n_medium=1.33, n_coverslip=1.33, n_immersion=1.33)

    @testset "Grazing collection requires a physical aperture" begin
        for args in ((), fixed), NA in (1.33, nextfloat(1.33), 1.4)
            error = try
                VectorPSF(NA, 0.68, args...; optics..., grid_size=33)
            catch caught
                caught
            end
            @test error isa ArgumentError
            if error isa ArgumentError
                @test occursin("90 degrees", error.msg)
                @test occursin("real objectives", error.msg)
                @test occursin("0.995 * n_medium", error.msg)
            end
        end
        vp = VectorPupilFunction(1.33, 0.68, 1.33, 1.33, 1.52, 33)
        @test_throws ArgumentError M.fill_vector_pupils!(vp, only(fixed))
        @test_throws ArgumentError M.fill_vector_pupils!(vp, only(fixed); normalize=false)
        for args in ((), fixed)
            for N in (32,33,64,65,128,129,255,256)
                limit = 0.995 * optics.n_medium
                a = VectorPSF(limit, 0.68, args...; optics..., grid_size=N)
                b = VectorPSF(prevfloat(limit), 0.68, args...;
                    optics..., grid_size=N)
                @test a(0.1, -0.2, 0.3) ≈ b(0.1, -0.2, 0.3) rtol=1e-10
                @test review_vector_energy(a) ≈ 1
                @test review_vector_energy(b) ≈ 1
                @test_throws ArgumentError VectorPSF(nextfloat(limit), 0.68, args...;
                    optics..., grid_size=N)
            end
        end
    end

    @testset "Finite grazing Fresnel-apodization product" begin
        # An axial rim sample has cm=0, while cc=ci>0. Its analytic coefficients
        # are Tp*A=2/sqrt(cc) and Ts*A=2*nm/(nc*sqrt(cc)).
        nm, nc, λ = 1.0, 1.5, 0.5
        d = DipoleVector(0.0, 0.0, 1.0)
        vp = VectorPupilFunction(nm, λ, nm, nc, nc, 5)
        M.fill_vector_pupils!(vp, d; normalize=false)
        cc = sqrt(1 - nm^2/nc^2)
        expected = 2 / sqrt(cc)
        @test vp.Ex.field[3,5] ≈ expected
        @test vp.Ey.field[3,5] ≈ 0 atol=1e-14
        transverse = VectorPupilFunction(nm, λ, nm, nc, nc, 5)
        M.fill_vector_pupils!(transverse, DipoleVector(0., 1., 0.); normalize=false)
        @test transverse.Ey.field[3,5] ≈ 2*nm/(nc*sqrt(cc))
        # The same critical angle can be strictly inside a supercritical aperture.
        inside = VectorPupilFunction(1.5, λ, nm, nc, nc, 7)
        M.fill_vector_pupils!(inside, d; normalize=false)
        @test inside.Ex.field[4,6] ≈ expected
        for NA in (prevfloat(nm), nm, nextfloat(nm))
            p = VectorPSF(NA, λ, d; n_medium=nm, n_coverslip=nc,
                n_immersion=nc, grid_size=5)
            @test isfinite(p(0.1, -0.2, 0.3))
            @test review_vector_energy(p) ≈ 1
        end
        p0 = VectorPSF(nm, λ, d; n_medium=nm, n_coverslip=nc,
            n_immersion=nc, grid_size=5)
        limit = p0(0.1, -0.2, 0.3)
        for sign in (-1, 1)
            errors = Float64[]
            for h in (1e-8, 1e-10, 1e-12)
                p = VectorPSF(nm*(1+sign*h), λ, d; n_medium=nm, n_coverslip=nc,
                    n_immersion=nc, grid_size=5)
                push!(errors, abs(p(0.1, -0.2, 0.3)-limit))
            end
            # The medium cosine is O(sqrt(h)); check approach to the finite limit.
            @test all(diff(errors) .< 0)
            @test errors[end] < 1e-4*abs(limit)
        end
    end

    @testset "Supercritical continuity and unit energy" begin
        for args in ((), fixed)
            values = Float64[]
            for N in (32, 64, 128)
                p = VectorPSF(1.4, 0.68, args...; n_medium=1.33,
                    n_coverslip=1.52, n_immersion=1.52, grid_size=N)
                @test review_vector_energy(p) ≈ 1 atol=2e-13
                value = p(0.1, -0.2, 0.3)
                @test isfinite(value)
                push!(values, value)
                for NA in (1.4-1e-8, 1.4+1e-8)
                    near = VectorPSF(NA, 0.68, args...; n_medium=1.33,
                        n_coverslip=1.52, n_immersion=1.52, grid_size=N)
                    @test review_vector_energy(near) ≈ 1 atol=2e-13
                    @test near(0.1, -0.2, 0.3) ≈ value rtol=1e-5
                end
            end
            @test maximum(abs.(diff(values))) / maximum(values) < 0.1
            @test abs(values[3]-values[2]) < abs(values[2]-values[1])
        end
    end

    @testset "Vector aperture agrees across all operations" begin
        for args in ((), fixed), N in (11, 33, 32, 64, 128)
            p = VectorPSF(1.2, 0.68, args...; optics..., grid_size=N)
            @test review_vector_energy(p) ≈ 1 atol=3e-13
            dk = 2*p.nₐ/p.λ/(N-1)
            for vp in p.vector_pupils
                shifted = copy(vp)
                M.apply_defocus!(shifted, 0.3)
                expected = ComplexF64[0, 0]
                outside_zero = outside_unchanged = true
                for i in 1:N, j in 1:N
                    kx, ky = (i-(N+1)/2)*dk, (j-(N+1)/2)*dk
                    kr2 = kx^2 + ky^2
                    inside = (2i - N - 1)^2 + (2j - N - 1)^2 <= (N - 1)^2    # exact aperture
                    if inside
                        kz = sqrt(complex((p.n_medium/p.λ)^2-kr2))
                        phase = cis(2π*(0.1*kx - 0.2*ky + 0.3*real(kz)))
                        expected[1] += vp.Ex.field[j,i]*phase*dk^2
                        expected[2] += vp.Ey.field[j,i]*phase*dk^2
                    else
                        outside_zero &= iszero(vp.Ex.field[j,i]) &&
                            iszero(vp.Ey.field[j,i])
                        outside_unchanged &= shifted.Ex.field[j,i] == vp.Ex.field[j,i] &&
                            shifted.Ey.field[j,i] == vp.Ey.field[j,i]
                    end
                end
                @test outside_zero
                @test outside_unchanged
                actual = M._calculate_field_amplitude(vp, p.nₐ, p.λ, p.n_medium,
                    p.n_coverslip, p.n_immersion, 0., 0.1, -0.2, 0.3)
                @test actual ≈ expected rtol=1e-12
                @test M._calculate_field_amplitude(shifted, p.nₐ, p.λ, p.n_medium,
                    p.n_coverslip, p.n_immersion, 0., 0.1, -0.2, 0.) ≈ expected rtol=1e-12
            end
        end
        # Ignored out-of-aperture storage must not reduce evaluated energy.
        p = VectorPSF(1.2, 0.68, only(fixed); optics..., grid_size=11)
        pupil = only(p.vector_pupils)
        pupil.Ex.field[1,1] = 1e6
        M.normalize!(pupil)
        @test review_vector_energy(p) ≈ 1 atol=2e-13
    end

    @testset "Radial polynomial integer widths" begin
        for I in (Int8, Int16, UInt8), (n,m) in ((16,12), (40,36), (40,0), (21,1))
            for r in (0.0, 0.5, 0.9, 1.0)
                @test radialpolynomial(I(n), I(m), r) == radialpolynomial(n, m, r)
            end
        end
        @test radialpolynomial(Int8(16), Int8(12), 0.5) == 0.01123046875
        for I in (Int16, UInt8), r in (0.5, 1.0)
            @test radialpolynomial(I(200), I(196), r) == radialpolynomial(200, 196, r)
        end
    end

    @testset "Nonfinite dipoles are rejected on every path" begin
        for bad in (Inf, -Inf, NaN), axis in 1:3
            c = zeros(3)
            c[axis] = bad
            @test_throws ArgumentError DipoleVector(c...)
            @test_throws ArgumentError DipoleVector{Float64}(c...)
            @test_throws ArgumentError DipoleVector{Float32}(c...)
            @test_throws ArgumentError DipoleVector{Int}(c...)
            @test_throws ArgumentError DipoleEmitter3D(1., 1., 0., 1., c...)
        end
        @test_throws ArgumentError DipoleVector(1., Float32(Inf), 0)
        @test DipoleVector{Float32}(floatmax(Float64), 0, 0).px === 1f0
        z = DipoleVector(0., 0., 0.)
        @test (z.px, z.py, z.pz) == (0., 0., 0.)
        @test DipoleVector{Int}(1, 0, 0) isa DipoleVector{Int}
    end

    @testset "Pupil coefficient construction keeps the evaluator matrix" begin
        c = ZernikeCoefficients(15)
        c.mag[6] = 0.2
        c.phase[4] = 0.3
        evaluate() = M.Zernike._evaluate_pupil(c, 128, Float64)
        construct() = PupilFunction(1.2, 0.68, 1.33, c; grid_size=128)
        @test construct().field == evaluate()
        evaluate(); construct()
        eval_bytes = minimum([@allocated(evaluate()) for _ in 1:3])
        make_bytes = minimum([@allocated(construct()) for _ in 1:3])
        @test make_bytes <= eval_bytes + 4096
    end
end
