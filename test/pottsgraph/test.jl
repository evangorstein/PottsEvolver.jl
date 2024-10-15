@testset "Gauge change" begin
    g = PottsGraph(5, 3; init=:rand)
    PottsEvolver.set_gauge!(g, :zero_sum)
    for i in 1:5
        @test abs(sum(g.h[:, i])) < 1e-5
        for j in 1:5, b in 1:3
            @test abs(sum(g.J[:, b, i, j])) < 1e-5
        end
    end

    @test_throws ErrorException PottsEvolver.set_gauge!(g, :some_random_gauge)
end
