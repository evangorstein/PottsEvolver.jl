using BioSequenceMappings

@testset "SamplingParameters" begin
    @test SamplingParameters(Teq=5) isa Any
    @test_throws AssertionError SamplingParameters(; Teq=5, step_type = :dubstep)
    @test_throws AssertionError SamplingParameters(; Teq=5, branch_length_type = :olive_tree)
end

@testset "Accepted steps" begin
    L, q, M = (5, 3, 100)
    g = PottsGraph(L, q; init = :rand)
    S, _ = mcmc_sample(g, M; branch_length_type = :accepted, Teq=1)
    # the difference between two samples should always be exactly one.
    @test map(i -> hamming(S[i-1], S[i]; normalize=false), 2:M) |> unique |> first == 1
end
