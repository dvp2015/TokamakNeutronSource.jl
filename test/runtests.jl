using TokamakNeutronSource
using Test

split(x, δ = 0.1) = [x - δ, x + δ]

@testset "TokamakNeutronSource.jl" begin
    # include("reaction-rates.jl")
    @testset "PlasmaDistribution" begin
        let excel_path = joinpath(@__DIR__, "data", "TRT_215_8T_NBI.xlsx")
            eqdsk_path = joinpath(@__DIR__, "data", "beforeTQ.eqdsk")
            @test isfile(eqdsk_path)
            using TokamakNeutronSource.PlasmaDistribution
            df = load_excel(excel_path)
            @testset "Excel" begin
                @test df[1, 1] == 0.0
                @test df[end, 1] == 1.0
            end
            eqdsk = Content(eqdsk_path)
            @testset "DD Distribution" begin
                distr = DDDistribution(eqdsk, df)
                @test distr.n(0) ≈ 9.436e13
                @test distr.n(2) == 0.0
                @test all(distr.n([0, 2]) .≈ [9.436e13, 0.0])
                @test Ti(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 21.09
                @test n(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 9.436e13
                @test I(distr, 0) ≈ 1.278e10 rtol=0.001
                @test I(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 1.278e10 rtol=0.001
                actual = I(distr, split(eqdsk.rmaxis), split(eqdsk.zmaxis))
                @test size(actual) == (2, 2)
            end
            @testset "DT Distribution" begin
                distr = DTDistribution(eqdsk, df)
                nd, nt = concentrations(distr, 0)
                @assert nd == nt == 0.5 * 9.436e13
                @test I(distr, 0) ≈ 1.034e12 rtol=0.001
                @test I(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 1.034e12 rtol=0.001
                # test vectorization
                actual = I(distr, split(eqdsk.rmaxis), split(eqdsk.zmaxis))
                @test size(actual) == (2, 2)
            end
        end
    end
end
