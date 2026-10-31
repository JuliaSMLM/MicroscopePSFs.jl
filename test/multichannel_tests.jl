using Random
using Statistics

# Probe comparators use only the public field interface.
struct ProbeScalarFourPi{A,B}
    top::A
    bottom::B
    phase::Float64
end
(p::ProbeScalarFourPi)(x, y, z) = abs2((amplitude(p.top, x, y, z) + cis(p.phase) * amplitude(p.bottom, x, y, -z)) / sqrt(2))
struct PolarizedPSF{P}
    psf::P
    component::Int
end
(p::PolarizedPSF)(x, y, z) = sum(abs2(E[p.component]) for E in dipole_fields(p.psf, x, y, z))

# Independent lab-frame oracle, ported from t27-mirror-check.jl. No package physics.
reflect(n) = (u = normalize(n); Matrix{Float64}(I, 3, 3) - 2 * u * u')
function rotate_onto(a, b)
    v = cross(a, b)
    c, s² = dot(a, b), dot(v, v)
    s² < 1e-24 && return Matrix{Float64}(I, 3, 3)
    K = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
    return Matrix{Float64}(I, 3, 3) + K + K * K * (1 - c) / s²
end
function arm_field(p, ze, xy; NA, λ, n, sgn, Nθ=240, Nϕ=240)
    k = 2π * n / λ
    θmax = asin(NA / n)
    E = zeros(ComplexF64, 3)
    dθ, dϕ = θmax / Nθ, 2π / Nϕ
    r = [-xy[1], -xy[2], ze]
    for iθ in 1:Nθ, iϕ in 1:Nϕ
        θ, ϕ = (iθ - 0.5) * dθ, (iϕ - 0.5) * dϕ
        k̂ = [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), sgn * cos(θ)]
        Efar = p - dot(p, k̂) * k̂
        Epup = rotate_onto(k̂, [0.0, 0.0, sgn]) * Efar
        E .+= Epup .* exp(-im * k * dot(k̂, r)) .* (sin(θ) * sqrt(cos(θ)))
    end
    return E .* dθ .* dϕ
end
function shape_err(code_vals, lab_vals)
    s = dot(code_vals, lab_vals) / dot(lab_vals, lab_vals)
    return maximum(abs.(code_vals .- s .* lab_vals)) / maximum(code_vals)
end

@testset "Multi-channel PSFs" begin
    rng = MersenneTwister(27)
    NA, λ, n = 1.2, 0.68, 1.33
    vector(d...; grid_size=256) = VectorPSF(NA, λ, d...; n_medium=n, n_coverslip=n, n_immersion=n, grid_size)
    fixed = vector(DipoleVector(1.0, 0.0, 0.0))
    rotating = vector()
    scalar = ScalarPSF(NA, λ, n)
    pts = [(0.0, 0.0, 0.1), (0.15, 0.05, 0.2), (-0.1, 0.2, -0.25), (0.25, -0.1, 0.35)]
    random_pts = [Tuple(rand(rng, 3) .* 0.4 .- 0.2) for _ in 1:8]
    xr, zr = -0.4:0.1:0.4, -0.4:0.1:0.4
    spline = SplinePSF(scalar, xr, xr, zr)
    spline2d = SplinePSF(GaussianPSF(0.15), xr, xr)

    @testset "dipole_fields" begin
        for pos in pts
            @test only(dipole_fields(fixed, pos...)) == amplitude(fixed, pos...)
            for p in (fixed, rotating, scalar)
                @test p(pos...) ≈ sum(sum(abs2, E) for E in dipole_fields(p, pos...)) rtol=1e-12
            end
            @test_throws ArgumentError dipole_fields(spline, pos...)
        end
        @test length(dipole_fields(rotating, pts[1]...)) == 3
        @test only(only(dipole_fields(scalar, pts[1]...))) == amplitude(scalar, pts[1]...)
    end

    @testset "AnalyzerPSF" begin
        for p in (fixed, rotating)
            h, v = AnalyzerPSF(p; component=1), AnalyzerPSF(p; component=2)
            for pos in pts
                @test h(pos...) + v(pos...) ≈ p(pos...) rtol=1e-12
                for c in 1:2
                    a = AnalyzerPSF(p; component=c)
                    @test a(pos...) ≈ PolarizedPSF(p, c)(pos...) rtol=1e-12
                    @test all(E -> iszero(E[3-c]), dipole_fields(a, pos...))
                    @test a(pos...) ≈ sum(sum(abs2, E) for E in dipole_fields(a, pos...)) rtol=1e-12
                end
            end
        end
        h = AnalyzerPSF(fixed; component=1)
        @test amplitude(h, pts[2]...) == only(dipole_fields(h, pts[2]...))
        err = try amplitude(rotating, pts[1]...) catch e; e end
        aerr = try amplitude(AnalyzerPSF(rotating; component=1), pts[1]...) catch e; e end
        @test sprint(showerror, err) == sprint(showerror, aerr)
        @test_throws ArgumentError AnalyzerPSF(scalar; component=1)
        @test_throws ArgumentError AnalyzerPSF(spline; component=1)
        for c in (0, 3, -1)
            @test_throws ArgumentError AnalyzerPSF(fixed; component=c)
        end
        # Pixel moments pin the orthogonal H/V squeezing without relying on a plot.
        fast = vector(; grid_size=64)
        edges = -0.8:0.08:0.8
        camera, emitter = IdealCamera(collect(edges), collect(edges)), Emitter3D(0.0, 0.0, 0.0, 1.0)
        centers = collect(edges[1:end-1]) .+ step(edges)/2
        images = [integrate_pixels(AnalyzerPSF(fast; component=c), camera, emitter; sampling=1) for c in 1:2]
        moments(img) = (sum(img .* reshape(centers.^2, 1, :)) / sum(img), sum(img .* centers.^2) / sum(img))
        hx, hy = moments(images[1]); vx, vy = moments(images[2])
        @test (hx - hy) * (vx - vy) < 0
        @test abs(hx - hy) > 0.01 * max(hx, hy)
        @test hx ≈ vy rtol=1e-12
        @test hy ≈ vx rtol=1e-12
    end

    @testset "PlanePSF" begin
        for p in (scalar, fixed, rotating, spline), pos in random_pts
            plane = PlanePSF(p; dz=0.1)
            x, y, z = pos
            @test plane(pos...) == p(x, y, z+plane.dz)
            if p !== spline
                @test dipole_fields(plane, pos...) == dipole_fields(p, x, y, z+plane.dz)
                p !== rotating && @test amplitude(plane, pos...) == amplitude(p, x, y, z+plane.dz)
            end
        end
        @test_throws ArgumentError PlanePSF(GaussianPSF(0.15); dz=0.1)
        @test_throws ArgumentError PlanePSF(spline2d; dz=0.1)
        @test_throws ArgumentError PlanePSF(scalar; dz=Inf)
    end

    @testset "Wrapper amplitude integration preserves polarization" begin
        fast = vector(DipoleVector(1.0, 0.0, 0.0); grid_size=64)
        camera = IdealCamera(4, 3, 0.1)
        emitter = Emitter3D(0.18, 0.14, 0.03, 1.0)
        emitters = [emitter, Emitter3D(0.22, 0.12, -0.02, 4.0)]
        plane = PlanePSF(fast; dz=0.1)
        analyzer = AnalyzerPSF(fast; component=1)
        cavity = FourPiPSF(fast, fast; phase=0.3)
        nested = PlanePSF(AnalyzerPSF(cavity; component=2); dz=0.05)
        kwargs = (; sampling=2, threaded=false)
        base = integrate_pixels_amplitude(fast, camera, emitter; kwargs...)
        shifted = Emitter3D(emitter.x, emitter.y, emitter.z+plane.dz, emitter.photons)
        @test integrate_pixels_amplitude(plane, camera, emitter; kwargs...) ==
              integrate_pixels_amplitude(fast, camera, shifted; kwargs...)
        analyzed = integrate_pixels_amplitude(analyzer, camera, emitter; kwargs...)
        @test analyzed[:,:,1] == base[:,:,1]
        @test all(iszero, analyzed[:,:,2])
        mirrored = Emitter3D(emitter.x, emitter.y, -emitter.z, emitter.photons)
        @test integrate_pixels_amplitude(cavity, camera, emitter; kwargs...) ≈
              (base .+ cis(cavity.phase) .* integrate_pixels_amplitude(fast, camera, mirrored; kwargs...)) ./ sqrt(2) rtol=1e-12
        for p in (plane, analyzer, cavity, nested)
            @test size(integrate_pixels_amplitude(p, camera, emitter; kwargs...)) == (3,4,2)
            combined = integrate_pixels_amplitude(p, camera, emitters; kwargs...)
            @test size(combined) == (3,4,2)
            @test combined ≈ sum(integrate_pixels_amplitude(p, camera, e; kwargs...) for e in emitters) rtol=1e-12
            empty_image = integrate_pixels_amplitude(p, camera, empty(emitters); kwargs...)
            @test size(empty_image) == (3,4,2)
            @test all(iszero, empty_image)
            for support in (0.2, (-0.1,0.1,-0.1,0.1)), source in (emitter, emitters)
                @test_throws ArgumentError integrate_pixels_amplitude(p, camera, source; support, kwargs...)
            end
        end
        for p in (PlanePSF(rotating; dz=0.1), AnalyzerPSF(rotating; component=1), FourPiPSF(rotating,rotating)),
            source in (emitter, emitters, empty(emitters))
            base_error = try integrate_pixels_amplitude(rotating, camera, source; kwargs...) catch e; e end
            wrapper_error = try integrate_pixels_amplitude(p, camera, source; kwargs...) catch e; e end
            @test wrapper_error isa ErrorException
            @test sprint(showerror, wrapper_error) == sprint(showerror, base_error)
        end
        for p in (PlanePSF(scalar; dz=0.1), FourPiPSF(scalar,scalar; phase=0.3)),
            source in (emitter, emitters, empty(emitters)),
            support in (Inf, 0.05, (0.1,0.2,0.1,0.2)), threaded in (false, true)
            actual = integrate_pixels_amplitude(p, camera, source; support, sampling=1, threaded)
            expected = integrate_pixels_amplitude(p, camera.pixel_edges_x, camera.pixel_edges_y,
                                                 source; support, sampling=1, threaded)
            @test actual isa Matrix
            @test size(actual) == (3,4)
            @test actual == expected
        end
        ambiguities = Test.detect_ambiguities(MicroscopePSFs)
        @test isempty(filter(pair -> any(m -> m.name == :integrate_pixels_amplitude, pair), ambiguities))
    end

    @testset "FourPiPSF scalar" begin
        ca, cb = ZernikeCoefficients(15), ZernikeCoefficients(15)
        ca.phase[6], cb.phase[8] = 0.3, 0.2
        top, bottom = ScalarPSF(NA, λ, n; zernike_coeffs=ca), ScalarPSF(NA, λ, n; zernike_coeffs=cb)
        for pos in random_pts
            x, y, z = pos
            @test FourPiPSF(top, bottom)(pos...) + FourPiPSF(top, bottom; phase=π)(pos...) ≈ top(pos...) + bottom(x, y, -z) rtol=1e-12
            for φ in (0.0, 0.3, π)
                @test FourPiPSF(top, bottom; phase=φ)(pos...) ≈ ProbeScalarFourPi(top, bottom, Float64(φ))(pos...) rtol=1e-12
            end
            p = FourPiPSF(PlanePSF(top; dz=0.1), bottom; phase=0.3)
            E = (amplitude(top, x, y, z+0.1) + cis(0.3)*amplitude(bottom, x, y, -z))/sqrt(2)
            @test amplitude(p, pos...) ≈ E rtol=1e-12
            @test p(pos...) ≈ abs2(E) rtol=1e-12
        end
        # Axial integral from probe2.jl, including its high-NA cos(theta) weight.
        naₐ, λₐ, nₐ = 1.35, 0.68, 1.516
        arm = ScalarPSF(naₐ, λₐ, nₐ)
        fp = FourPiPSF(arm, arm)
        zs = collect(-0.6:0.004:0.6)
        us = range(sqrt(1-(naₐ/nₐ)^2), 1; length=4001)
        analytic = [abs2(sum(cos(2π*(nₐ/λₐ)*u*z)*u for u in us)) for z in zs]
        curve = [fp(0.0, 0.0, z) for z in zs]
        curve ./= maximum(curve); analytic ./= maximum(analytic)
        peaks(a) = [i for i in 2:length(a)-1 if a[i] > a[i-1] && a[i] > a[i+1]]
        @test abs(median(diff(zs[peaks(curve)])) - median(diff(zs[peaks(analytic)]))) <= 0.004 + eps()
        @test maximum(abs.(curve .- analytic)) < 0.01
    end

    @testset "FourPiPSF vector, independent cavity oracle" begin
        Rtop = reflect([1, 0, 1])
        Rbot = reflect([1, 0, 0]) * reflect([1, 0, 1])
        @test Rtop[:, 1:2] ≈ Rbot[:, 1:2]
        @test Rtop * [0, 0, -1] ≈ Rbot * [0, 0, 1]
        lab_dipoles = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0], normalize([1.0, 0, -1.0])]
        # Compute each arm only once per dipole and position; reuse across cavity phases.
        lab_arms = [[(Rtop * arm_field(p, z, (x,y); NA, λ, n, sgn=-1),
                      Rbot * arm_field(p, z, (x,y); NA, λ, n, sgn=1)) for (x,y,z) in pts] for p in lab_dipoles]
        code_dipoles = [DipoleVector(1.0,0.0,0.0), DipoleVector(0.0,1.0,0.0), DipoleVector(0.0,0.0,1.0), DipoleVector(1.0,0,1.0)]
        code_pairs = [(vector(d), vector(DipoleVector(d.px, d.py, -d.pz))) for d in code_dipoles]
        for (idx, (top, bottom)) in enumerate(code_pairs), φ in (0.0, π/2, π)
            lab_vals = [sum(abs2, a .+ cis(φ).*b) / 2 for (a,b) in lab_arms[idx]]
            p = FourPiPSF(top, bottom; phase=φ)
            @test shape_err([p(pos...) for pos in pts], lab_vals) < 0.002
            @test amplitude(p, pts[2]...) == only(dipole_fields(p, pts[2]...))
        end
        for φ in (0.0, π/2, π)
            lab_vals = [sum(sum(abs2, lab_arms[d][j][1] .+ cis(φ).*lab_arms[d][j][2])/2 for d in 1:3) for j in eachindex(pts)]
            p = FourPiPSF(rotating, rotating; phase=φ)
            @test shape_err([p(pos...) for pos in pts], lab_vals) < 0.002
            for (x,y,z) in pts
                @test p(x,y,z) + FourPiPSF(rotating, rotating; phase=φ+π)(x,y,z) ≈ rotating(x,y,z) + rotating(x,y,-z) rtol=1e-12
            end
        end
        # Omitting the axial-dipole mirror sign must fail, independently of the implementation.
        no_sign = [sum(abs2, dipole_fields(rotating,x,y,z)[3] .+ dipole_fields(rotating,x,y,-z)[3])/2 for (x,y,z) in pts]
        lab_z = [sum(abs2, a .+ b)/2 for (a,b) in lab_arms[3]]
        @test shape_err(no_sign, lab_z) > 0.5
        # Existing VectorPSF frame: pz points opposite emitter z. Code (1,0,1)
        # therefore matches lab (1,0,-1), not lab (1,0,1). Mirroring pz holds in either frame.
        tilt = code_pairs[4][1]
        code_vals = [tilt(pos...) for pos in pts]
        @test shape_err(code_vals, [sum(abs2,a) for (a,b) in lab_arms[4]]) < 0.002
        wrong_lab = [sum(abs2,arm_field(normalize([1.0,0,1.0]),z,(x,y); NA,λ,n,sgn=-1)) for (x,y,z) in pts]
        @test shape_err(code_vals, wrong_lab) > 0.1
        err = try amplitude(rotating, pts[1]...) catch e; e end
        ferr = try amplitude(FourPiPSF(rotating,rotating), pts[1]...) catch e; e end
        @test sprint(showerror, err) == sprint(showerror, ferr)
    end

    @testset "4Pi refusals" begin
        evanescent = VectorPSF(1.4, λ; n_medium=n, grid_size=64)
        @test_throws ArgumentError FourPiPSF(evanescent, evanescent)
        tilt = vector(DipoleVector(1.0,0.0,1.0))
        axial = vector(DipoleVector(0.0,0.0,1.0))
        @test_throws ArgumentError FourPiPSF(tilt, tilt)
        @test_throws ArgumentError FourPiPSF(axial, axial)
        @test_throws ArgumentError FourPiPSF(PlanePSF(tilt; dz=0.1), AnalyzerPSF(tilt; component=1))
        @test_throws ArgumentError FourPiPSF(scalar, fixed)
        @test_throws ArgumentError FourPiPSF(fixed, rotating)
        @test_throws ArgumentError FourPiPSF(spline, spline)
        @test_throws ArgumentError FourPiPSF(scalar, scalar; phase=NaN)
        # Normalize before comparison even for the same-typed constructor path.
        a = vector(DipoleVector(1.0,1.0,0.0))
        b = vector(DipoleVector(2.0,2.0,0.0))
        @test FourPiPSF(a,b) isa FourPiPSF
        mirror = vector(DipoleVector(1.0,0.0,-1.0))
        @test FourPiPSF(AnalyzerPSF(PlanePSF(tilt; dz=0.1); component=1), AnalyzerPSF(mirror; component=1)) isa FourPiPSF
    end

    @testset "4Pi rejects malformed pupil counts" begin
        malformed = vector(; grid_size=64)
        pop!(malformed.vector_pupils)
        @test length(malformed.vector_pupils) == 2
        @test_throws ArgumentError FourPiPSF(malformed, malformed)
        err = try FourPiPSF(malformed, malformed) catch e; e end
        @test err isa ArgumentError
        if err isa ArgumentError
            @test occursin("2 pupils", err.msg)
        end
        empty!(malformed.vector_pupils)
        @test_throws ArgumentError FourPiPSF(malformed, malformed)
        append!(malformed.vector_pupils, fill(first(rotating.vector_pupils), 4))
        @test_throws ArgumentError FourPiPSF(malformed, malformed)
    end

    @testset "4Pi requires one optical frequency" begin
        other = ScalarPSF(NA, 0.80, n)
        @test_throws ArgumentError FourPiPSF(scalar, other)
        @test_throws ArgumentError FourPiPSF(PlanePSF(scalar; dz=0.1), PlanePSF(other; dz=-0.1))
        other_vector = VectorPSF(NA, 0.80, DipoleVector(1.0,0.0,0.0);
                                 n_medium=n, n_coverslip=n, n_immersion=n, grid_size=64)
        @test_throws ArgumentError FourPiPSF(AnalyzerPSF(fixed; component=1),
                                            PlanePSF(AnalyzerPSF(other_vector; component=1); dz=0.1))
        err = try FourPiPSF(scalar, other) catch e; e end
        @test err isa ArgumentError
        if err isa ArgumentError
            @test occursin(string(scalar.λ), err.msg)
            @test occursin(string(other.λ), err.msg)
        end
        τ = sqrt(eps(Float64))
        @test FourPiPSF(scalar, ScalarPSF(NA, λ*(1+τ/2), n)) isa FourPiPSF
        @test_throws ArgumentError FourPiPSF(scalar, ScalarPSF(NA, λ*(1+2τ), n))
        @test FourPiPSF(scalar, ScalarPSF(1.1, λ, n)) isa FourPiPSF
    end

    @testset "4Pi tolerates mixed arm precision" begin
        scalar32 = ScalarPSF(1.2f0, 0.66f0, 1.33f0; pupil_data=ones(ComplexF32,16,16))
        scalar64 = ScalarPSF(1.2, 0.66, 1.33; pupil_data=ones(ComplexF64,16,16))
        @test scalar32 isa ScalarPSF{Float32}
        @test FourPiPSF(scalar32, scalar64) isa FourPiPSF{Float64}
        @test FourPiPSF(scalar64, scalar32) isa FourPiPSF{Float64}
        @test FourPiPSF(PlanePSF(scalar32; dz=0.1), scalar64) isa FourPiPSF{Float64}
        @test_throws ArgumentError FourPiPSF(
            ScalarPSF(1.2f0, 0.68f0, 1.33f0; pupil_data=ones(ComplexF32,16,16)),
            ScalarPSF(1.2, 0.80, 1.33; pupil_data=ones(ComplexF64,16,16)))

        # Preserve physical fields in Float32 pupils; the vector factory allocates Float64 pupils.
        function float32_arm(p::VectorPSF)
            nₐ, λ, nₘ, n꜀, nᵢ = Float32.((p.nₐ,p.λ,p.n_medium,p.n_coverslip,p.n_immersion))
            pupils = [VectorPupilFunction(nₐ,λ,nₘ,n꜀,nᵢ,
                PupilFunction(nₐ,λ,nₘ,ComplexF32.(v.Ex.field)),
                PupilFunction(nₐ,λ,nₘ,ComplexF32.(v.Ey.field))) for v in p.vector_pupils]
            d = p.dipole
            return VectorPSF{Float32}(nₐ,λ,nₘ,n꜀,nᵢ,
                DipoleVector(Float32(d.px),Float32(d.py),Float32(d.pz)),
                Float32(p.z_stage),pupils,nothing,nothing)
        end
        # Exactly representable wavelength isolates the dipole tolerance from the wavelength check.
        arm(d) = VectorPSF(1.2, 0.5, d; n_medium=n,n_coverslip=n,n_immersion=n,grid_size=16)
        for (px,py,pz) in ((0.2,0.4,0.66), (1.0,2.0,3.0), (0.91,0.32,0.17))
            top32 = float32_arm(arm(DipoleVector(px,py,pz)))
            bottom64 = arm(DipoleVector(px,py,-pz))
            @test FourPiPSF(top32, bottom64) isa FourPiPSF{Float64}
            @test FourPiPSF(bottom64, top32) isa FourPiPSF{Float64}
            @test FourPiPSF(PlanePSF(top32; dz=0.1), AnalyzerPSF(bottom64; component=1)) isa FourPiPSF{Float64}
            @test_throws ArgumentError FourPiPSF(top32, arm(DipoleVector(px,py,pz)))
        end
    end

    @testset "Analyzer over 4Pi" begin
        p = FourPiPSF(rotating, rotating; phase=0.4)
        h, v = AnalyzerPSF(p; component=1), AnalyzerPSF(p; component=2)
        for pos in pts
            @test h(pos...) + v(pos...) ≈ p(pos...) rtol=1e-12
        end
    end

    @testset "Spline tabulation of wrappers" begin
        fast = vector(; grid_size=64)
        pupil = PupilFunction(NA, λ, n, ZernikeCoefficients(1); grid_size=64)
        sc = ScalarPSF(NA, λ, n; pupil)
        xx, zz = -0.35:0.05:0.35, -0.18:0.03:0.18
        wrappers = (PlanePSF(sc; dz=0.1), PlanePSF(fast; dz=0.1),
                    AnalyzerPSF(fast; component=1), FourPiPSF(sc,sc; phase=0.3), FourPiPSF(fast,fast; phase=0.3))
        interior = [Tuple((rand(rng,3).-0.5) .* [0.44,0.44,0.2]) for _ in 1:12]
        for p in wrappers
            s = SplinePSF(p,xx,xx,zz)
            peak = maximum(s.original_grid)
            @test maximum(abs(s(pos...) - p(pos...)) for pos in interior) < 0.02*peak
        end
    end

    @testset "MultiChannelPSF" begin
        mc = MultiChannelPSF(; channels=Any[spline,spline])
        @test !(mc isa AbstractPSF)
        @test eltype(mc.channels) == typeof(spline)
        @test mc.transmission == [0.5,0.5]
        @test mc.z_offset == [0.0,0.0]
        @test mc.reference == 1
        @test mc.registration == [[1.0 0 0; 0 1 0], [1.0 0 0; 0 1 0]]
        @test MultiChannelPSF(;channels=[scalar]).transmission == [1.0]
        @test MultiChannelPSF(;channels=[scalar],transmission=Float32[1]).transmission isa Vector{Float32}
        @test MultiChannelPSF(;channels=[scalar],transmission=Float32[1],z_offset=[0.0]).transmission isa Vector{Float64}
        @test_throws MethodError integrate_pixels(mc, IdealCamera(collect(0:0.1:1),collect(0:0.1:1)), Emitter3D(0.5,0.5,0.0,1.0))
        for kwargs in ((;channels=[]), (;channels=[scalar,1]),
                (;channels=[scalar],registration=[]), (;channels=[scalar],transmission=[]), (;channels=[scalar],z_offset=[]),
                (;channels=[scalar],registration=[ones(2,2)]), (;channels=[scalar],registration=[zeros(3,2)]),
                (;channels=[scalar],registration=[[1.0 0 Inf; 0 1 0]]),
                (;channels=[scalar],transmission=[NaN]), (;channels=[scalar],z_offset=[Inf]),
                (;channels=[scalar],transmission=[-0.1]), (;channels=[scalar,scalar],transmission=[1.0,1.0]),
                (;channels=[scalar],reference=0), (;channels=[scalar],reference=2), (;channels=[scalar],reference=1.5),
                (;channels=[scalar],registration=[[1.0 0 0.1;0 1 0]]), (;channels=[scalar],z_offset=[0.1]))
            @test_throws ArgumentError MultiChannelPSF(;kwargs...)
        end
        τ = sqrt(eps(Float64))
        @test MultiChannelPSF(;channels=[scalar],transmission=[1+τ/2]) isa MultiChannelPSF
        @test_throws ArgumentError MultiChannelPSF(;channels=[scalar],transmission=[1+2τ])
        @test MultiChannelPSF(;channels=[scalar],registration=[[1.0 0 τ/2;0 1 0]],z_offset=[τ/2]) isa MultiChannelPSF
        r = [0.0 -1 2; 1 0 3]
        registered = MultiChannelPSF(;channels=[scalar,scalar],registration=[mc.registration[1],r],z_offset=[0.0,0.2])
        @test collect(channel_position(registered,2,0.3,0.4,0.1)) ≈ [1.6,3.3,0.3]
        @test_throws ArgumentError channel_position(registered,3,0,0,0)
        @test MultiChannelPSF(;channels=[scalar,scalar],reference=2).reference == 2
        for p in (PlanePSF(scalar;dz=0.1), AnalyzerPSF(fixed;component=1), FourPiPSF(scalar,scalar), mc)
            @test occursin(string(nameof(typeof(p))), sprint(show,p))
        end

        # H/V tabulated as photon-per-grid-cell arrays, each scaled to unit sum.
        # Divide by cell area to recover a unit-integral intensity density.
        fast = vector(;grid_size=64)
        xx = -1.0:0.05:1.0
        area = step(xx)^2
        grids = [[AnalyzerPSF(fast;component=c)(x,y,0.0)*area for y in xx,x in xx] for c in 1:2]
        sums = sum.(grids)
        normalized = [SplinePSF((g ./ sum(g)) ./ area,xx,xx) for g in grids]
        hv = MultiChannelPSF(;channels=normalized, transmission=sums./sum(sums))
        cam, em = IdealCamera(collect(-0.7:0.1:0.7),collect(-0.7:0.1:0.7)), Emitter3D(0.0,0.0,0.0,1.0)
        expected = integrate_pixels(fast,cam,em; sampling=1) ./ sum(sums)
        N = 1200.0
        images = [N*hv.transmission[i]*integrate_pixels(hv.channels[i],cam,
                  Emitter3D(channel_position(hv,i,em.x,em.y,em.z)...,1.0);sampling=1) for i in 1:2]
        @test maximum(abs.(sum(images) .- N.*expected)) < 0.02 * maximum(N.*expected)
    end

    @testset "Container numeric types are writable" begin
        for T in (Float16, BigFloat)
            for kwargs in ((;transmission=T[1]), (;z_offset=T[0]),
                           (;registration=[T[1 0 0;0 1 0]]))
                @test_throws ArgumentError MultiChannelPSF(;channels=[scalar], kwargs...)
            end
            err = try MultiChannelPSF(;channels=[scalar], transmission=T[1]) catch e; e end
            @test err isa ArgumentError
            if err isa ArgumentError
                @test occursin("Float32 or Float64", err.msg)
            end
        end
        @test MultiChannelPSF(;channels=[scalar], registration=[Float32[1 0 0;0 1 0]]) isa MultiChannelPSF{Float32}
        @test MultiChannelPSF(;channels=[scalar], transmission=[1], z_offset=[0],
                             registration=[[1 0 0;0 1 0]]) isa MultiChannelPSF{Float64}
        @test MultiChannelPSF(;channels=[scalar], transmission=[1//2], z_offset=[0//1]) isa MultiChannelPSF{Float64}
        @test MultiChannelPSF(;channels=[scalar], transmission=[big(1)]) isa MultiChannelPSF{Float64}
    end

    @testset "Container refuses channels carrying a photon split" begin
        h, v = AnalyzerPSF(fixed; component=1), AnalyzerPSF(fixed; component=2)
        @test_throws ArgumentError MultiChannelPSF(;channels=[h,v])
        for p in (h, PlanePSF(h; dz=0.1), FourPiPSF(h,h),
                  FourPiPSF(fixed,h), FourPiPSF(h,fixed),
                  AnalyzerPSF(FourPiPSF(fixed,fixed); component=1))
            @test_throws ArgumentError MultiChannelPSF(;channels=[p])
        end
        err = try MultiChannelPSF(;channels=[h,v]) catch e; e end
        @test err isa ArgumentError
        if err isa ArgumentError
            @test occursin("AnalyzerPSF carries the photon split itself", err.msg)
            @test occursin("SplinePSF scaled to unit integral", err.msg)
            @test occursin("transmission", err.msg)
        end
        # Unit-integral spline channels are covered by the H/V forward-law test above.
    end

    @testset "Nested I/O and import-only child" begin
        mktempdir() do dir
            file = joinpath(dir,"multichannel.h5")
            mc = MultiChannelPSF(;channels=Any[spline,spline], reference=2,
                registration=[[0.0 -1 0.2;1 0 0.3],[1.0 0 0;0 1 0]], transmission=[0.3,0.6], z_offset=[0.1,0.0])
            save_psf(file,mc;metadata=Dict("description"=>"T27", "test_id"=>27))
            loaded = load_psf(file)
            @test loaded.channels isa Vector{typeof(spline)}
            @test loaded.registration == mc.registration
            @test loaded.transmission == mc.transmission
            @test loaded.z_offset == mc.z_offset
            @test loaded.reference == mc.reference
            @test all(loaded.channels[i](pos...) == mc.channels[i](pos...) for i in 1:2,pos in pts)
            h5open(file,"r") do f
                @test read(attributes(f)["psf_type"]) == "MultiChannelPSF"
                @test read(attributes(f)["description"]) == "T27"
                @test read(attributes(f)["test_id"]) == "27"
                @test size(read(f["parameters/registration"])) == (2,3,2)
                @test read(attributes(f["channels/channel_1"])["psf_type"]) == "SplinePSF"
            end
            script = """
                import MicroscopePSFs
                mc = MicroscopePSFs.load_psf(ARGS[1])
                @assert isfinite(mc.channels[1](0.1,0.0,0.1))
                @assert all(id.name ∉ ("Reactant", "PSFLearning") for id in keys(Base.loaded_modules))
                print("true")
                """
            cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script $file`
            out = IOBuffer()
            @test success(pipeline(cmd;stdout=out))
            @test String(take!(out)) == "true"
            wrappers = [PlanePSF(scalar;dz=0.1), AnalyzerPSF(rotating;component=1),
                        FourPiPSF(scalar,scalar;phase=0.3),
                        AnalyzerPSF(FourPiPSF(rotating,rotating;phase=0.4);component=1),
                        PlanePSF(AnalyzerPSF(rotating;component=1);dz=0.1)]
            all_channels = [scalar,fixed,rotating,wrappers[1],wrappers[3],FourPiPSF(rotating,rotating;phase=0.4)]
            heterogeneous = MultiChannelPSF(;channels=all_channels)
            save_psf(file,heterogeneous)
            restored = load_psf(file)
            for (a,b) in zip(restored.channels,all_channels), pos in pts
                @test a(pos...) == b(pos...)
            end
            for wrapper in wrappers
                save_psf(file,wrapper)
                restored = load_psf(file)
                @test typeof(restored) == typeof(wrapper)
                @test all(restored(pos...) == wrapper(pos...) for pos in pts)
            end
            # Loading reconstructs wrappers and the container through their checked constructors.
            save_psf(file,FourPiPSF(scalar,scalar))
            h5open(file,"r+") do f; write(f["parameters/phase"], NaN); end
            @test_throws ErrorException load_psf(file)
            save_psf(file,mc)
            h5open(file,"r+") do f; write(f["parameters/transmission"], [1.0,1.0]); end
            @test_throws ErrorException load_psf(file)
        end
    end

    @testset "Additive exports" begin
        # Hard-coded origin/main exports, including the calculate_* physics primitives.
        original = [:AbstractCamera,:AbstractEmitter,:AbstractPSF,:AiryPSF,:DipoleEmitter3D,:DipoleVector,
            :Emitter2D,:Emitter3D,:GaussianPSF,:IdealCamera,:MicroscopePSFs,:PupilFunction,:ScalarPSF,
            :SplinePSF,:VectorPSF,:VectorPupilFunction,:Zernike,:ZernikeCoefficients,:amplitude,
            :calculate_apodization,:calculate_axial_phase,:calculate_dipole_field_components,
            :calculate_fresnel_coefficients,:calculate_interface_fresnel,:calculate_pupil_field,
            :calculate_wave_vectors,:get_pixel_centers,:get_pixel_size,:integrate_pixels,
            :integrate_pixels_amplitude,:is_axis_aligned,:load_psf,:pupil_at,:radialpolynomial,
            :save_psf,:zernikepolynomial]
        @test all(n -> n in names(MicroscopePSFs), original)
        @test all(n -> n in names(MicroscopePSFs), [:dipole_fields,:PlanePSF,:AnalyzerPSF,:FourPiPSF,:MultiChannelPSF,:channel_position])
    end
end
