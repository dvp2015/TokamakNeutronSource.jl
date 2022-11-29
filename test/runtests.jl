using TokamakNeutronSource
using Test

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
            distr = Distribution(eqdsk, df)
            @testset "Distribution" begin
                @test T(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 21.09
                @test n(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 9.436e13
                @test Idd(distr, 0) ≈ 1.278e10 rtol=0.001
                @test Idd(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 1.278e10 rtol=0.001
                @test Idt(distr, 0) ≈ 1.034e12 rtol=0.001
                @test Idt(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 1.034e12 rtol=0.001
            end
        end
    end
end
