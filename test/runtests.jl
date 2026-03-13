using TokamakNeutronSource
using Test

split(x::Float64; δ::Float64=0.1)::Vector{Float64} = [x - δ, x + δ]

@testset "TokamakNeutronSource.jl" begin
    # include("reaction-rates.jl")
    @testset "PlasmaDistributions" begin
        let excel_path = joinpath(@__DIR__, "data", "TRT_215_8T_NBI.xlsx")
            eqdsk_path = joinpath(@__DIR__, "data", "beforeTQ.eqdsk")
            @test isfile(eqdsk_path)
            using TokamakNeutronSource.PlasmaDistributions
            using TokamakNeutronSource.Integrations
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
                @test I(distr, 0) ≈ 1.278e10 rtol = 0.001
                @test I(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 1.278e10 rtol = 0.001
                actual = I(distr, split(eqdsk.rmaxis), split(eqdsk.zmaxis))
                @test size(actual) == (2, 2)
                total, err, neval, fail = total_yield(distr)
                # 1-st moment of distribution
                rmin, rmax, zmin, zmax = domain(distr)
                @assert fail == 0
                @test total ≈ 9.939e16 rtol = 1e-4
                @test err / total < 1e-4
                (r1, z1), (r1err, z1err), neval, fail = torroidal_segment_moment_1(
                    (r, z) -> I(distr, r, z), (rmin, rmax), (zmin, zmax)
                )
                @assert fail == 0
                @test r1 ≈ 2.187 rtol = 3e-4
                @test z1 ≈ 0.496 rtol = 5e-4
                @test r1err / r1 < 3e-4
                @test z1err / z1 < 3e-4
            end
            @testset "DT Distribution" begin
                distr = DTDistribution(eqdsk, df)
                nd, nt = concentrations(distr, 0)
                @assert nd == nt == 0.5 * 9.436e13
                @test I(distr, 0) ≈ 1.034e12 rtol = 0.001
                @test I(distr, eqdsk.rmaxis, eqdsk.zmaxis) ≈ 1.034e12 rtol = 0.001
                # test vectorization
                actual = I(distr, split(eqdsk.rmaxis), split(eqdsk.zmaxis))
                @test size(actual) == (2, 2)
                total, err, neval, fail = total_yield(distr)
                @assert fail == 0
                @test total ≈ 8.438e18 rtol = 1e-4
                @test err / total < 1e-4
            end
        end
    end
    @testset "Testing utils" begin
        using TokamakNeutronSource.Integrations
        using TokamakNeutronSource.Testing
        total, err, neval, fail = total_yield(TestDistribution())
        @assert fail == 0
        @test total ≈ 0.5 * 1e6 * 2π * (2.5^2 - 1.5^2) rtol = 1e-4
        @test err / total < 1e-4
    end
end
