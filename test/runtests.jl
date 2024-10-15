using Accessors
using StatsBase
using PottsEvolver
using Test

@testset "PottsEvolver.jl" begin
    @testset "Codons" begin
        include("codons/test.jl")
    end

    @testset "PottsGraph" begin
        include("pottsgraph/test.jl")
    end

    @testset "Sequence" begin
        include("sequences/test.jl")
    end

    @testset "Sampling" begin
        include("sampling/test.jl")
    end
end
