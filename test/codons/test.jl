@testset "GenCode: Symbol to Int to Symbol" begin
    for i in 1:length(codon_alphabet)
        aa_char = genetic_code(codon_alphabet[i])
        aa_int = genetic_code(i)
        @test (isnothing(aa_int) && aa_char == '*') || (aa_alphabet(aa_int) == aa_char)
    end
end

@testset "Genetic code" begin
    # This mainly tests self-consistency
    for aa in 1:length(aa_alphabet)
        @test aa in map(genetic_code, PottsEvolver.reverse_code(aa))
        @test PottsEvolver.reverse_code_rand(aa) |> genetic_code == aa

        aa_c = aa_alphabet(aa)
        @test aa_c in map(genetic_code, PottsEvolver.reverse_code(aa_c))
        @test PottsEvolver.reverse_code_rand(aa_c) |> genetic_code == aa_c
    end
end

@testset "Valid/Gap/Invalid codons" begin
    stop_codons = ["TAA", "TAG", "TGA"]
    for codon in symbols(codon_alphabet)
        @test in(prod(bases(codon)), stop_codons) == PottsEvolver.isstop(codon)
        # Codons made of only nucleotides are always valid
        @test PottsEvolver.isvalid(codon) == (
            all(in(['A', 'C', 'G', 'T']), bases(codon)) ||
            all(==('-'), bases(codon))
        )
    end
end


