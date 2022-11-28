@testset "Reaction Rates" begin
    using TokamakNeutronSource.ReactionRates

    epsilon = 0.001

    temperatures = [
        0.0,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        1.0,
        # 1.3, # invalid value in the Bosch&Hale table
        1.5,
        # 1.8, #
        2.0,
        2.5,
        3.0,
        4.0,
        5.0,
        6.0,
        8.0,
        10.0,
        12.0,
        15.0,
        20.0,
        30.0,
        40.0,
        50.0,
    ]

    expected_dt = [
        0.0,
        1.254e-26,
        7.292e-25,
        9.344e-24,
        5.697e-23,
        2.253e-22,
        6.740e-22,
        1.662e-21,
        6.857e-21,
        # 2.546e-20, #  invalid value
        6.923e-20,
        # 1.539e-19, #  invalid value
        2.977e-19,
        8.425e-19,
        1.867e-18,
        5.974e-18,
        1.366e-17,
        2.554e-17,
        6.222e-17,
        1.136e-16,
        1.747e-16,
        2.740e-16,
        4.330e-16,
        6.681e-16,
        7.998e-16,
        8.649e-16,
    ]

    expected_ddn = [
        0.0,
        4.482e-28,
        2.004e-26,
        2.168e-25,
        1.169e-24,
        4.200e-24,
        1.162e-23,
        2.681e-23,
        9.933e-23,
        # 3.319e-22, # invalid value
        8.284e-22,
        # 1.713e-21,
        3.110e-21,
        7.905e-21,
        1.602e-20,
        4.447e-20,
        9.128e-20,
        1.573e-19,
        3.457e-19,
        6.023e-19,
        9.175e-19,
        1.481e-18,
        2.603e-18,
        5.271e-18,
        8.235e-18,
        1.133e-17,
    ]

    expected = hcat(expected_dt, expected_ddn)

    @testset "Selecting from DT and DD" begin
        for (i, collided) in enumerate([DT, DDN])
            @testset "Scalar values, $collided" begin
                for (t, e) in zip(temperatures, expected[:, i])
                    f = σv(collided)
                    @test f(t) ≈ e rtol = epsilon
                end
            end

            @testset "Vector argument, $collided" begin
                f = σv(collided)
                actual = f.(temperatures)
                @test maximum(abs.(actual .- expected[:, i])) < epsilon
            end

            @testset "Matix argument, $collided" begin
                indices = [1 3; 5 7]
                temps = temperatures[indices]
                f = σv(collided)
                actual = f.(temps)   # alias should also work
                @test maximum(abs.(actual .- expected[indices, i])) < epsilon
            end
        end
    end
end
