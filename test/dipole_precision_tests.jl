using Test, MicroscopePSFs

@testset "Dipole normalization working precision" verbose=true begin
    @testset "Mixed integer and Float16 inputs retain a finite direction" begin
        d = DipoleVector{Float64}(65536,Float16(0),Float16(0))
        @test d isa DipoleVector{Float64}
        @test (d.px,d.py,d.pz) === (1.0,0.0,0.0)
    end

    @testset "Integer inputs use the requested BigFloat precision" begin
        setprecision(BigFloat,256) do
            d = DipoleVector{BigFloat}(2,3,6)
            @test d isa DipoleVector{BigFloat}
            @test sum(abs2,(d.px,d.py,d.pz)) == one(BigFloat)
        end
    end
end
