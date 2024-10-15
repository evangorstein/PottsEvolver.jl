#====================================#
############# AASequence #############
#====================================#

@testset "AASequence Tests" begin
    # Test constructing the object
    aa_seq = AASequence([1, 2, 3, 4, 5])
    @test aa_seq isa AASequence
    @test length(aa_seq) == 5

    # Test indexing and setting the index in a sequence
    @test aa_seq[2] == 2
    aa_seq[2] = 10
    @test aa_seq[2] == 10

    # Test copying the object
    aa_seq_copy = copy(aa_seq)
    @test aa_seq_copy == aa_seq
    @test pointer(aa_seq.seq) != pointer(aa_seq_copy.seq)
end

#=============================================#
################ CodonSequence ################
#=============================================#

@testset "CodonSequence" begin
    @testset "From aa" begin
        aas = collect(1:6) # Int
        codon_seq = CodonSequence(aas; source=:aa)

        @test codon_seq isa CodonSequence{Int}
        @test PottsEvolver.translate(codon_seq).seq == codon_seq.aaseq
    end

    @testset "From codons" begin
        # Create an example codon sequence with matching amino acids
        codons = Int8.([1, 2, 3, 6, 7])
        aas = map(genetic_code, codons)
        codon_seq = CodonSequence(codons)

        # Test constructing the object
        @test codon_seq isa CodonSequence{Int8}
        @test length(codon_seq) == 5

        # Test indexing and setting the index in a sequence
        @test codon_seq[3] == 3
        codon_seq[3] = 7
        @test codon_seq[3] == 7
        @test codon_seq.aaseq[3] == genetic_code(7)

        # Test copying the object
        codon_seq_copy = copy(codon_seq)
        @test codon_seq_copy == codon_seq
        @test pointer(codon_seq.seq) != pointer(codon_seq_copy.seq)
        @test pointer(codon_seq.aaseq) != pointer(codon_seq_copy.aaseq)
    end
end

#============================================================#
##################### Numerical sequence #####################
#============================================================#

@testset "NumSequence Tests" begin
    # Test constructing the object
    L, q = (10, 5)
    num_seq = PottsEvolver.NumSequence(L, q; T=Int32)
    @test num_seq isa PottsEvolver.NumSequence{Int32}
    @test length(num_seq) == L

    # Test indexing and setting the index in a sequence
    @test num_seq[4] >= 1 && num_seq[4] <= q
    num_seq[4] = 2
    @test num_seq[4] == 2

    # Test copying the object
    num_seq_copy = copy(num_seq)
    @test num_seq_copy == num_seq
    @test pointer(num_seq.seq) != pointer(num_seq_copy.seq)
end

@testset "Sequences to alignent" begin
    M, L, q = (3, 10, 5)
    numseqs = map(_ -> PottsEvolver.NumSequence(L, q; T=Int32), 1:M)
    aaseqs = map(_ -> AASequence(L), 1:M) # default integer type should be IntType
    codonseqs = map(_ -> CodonSequence(L), 1:M)

    A = Alignment(numseqs)
    @test A isa Alignment{Nothing,Int32}
    @test isnothing(A.alphabet)

    A = Alignment(aaseqs)
    @test A isa Alignment{Char,PottsEvolver.IntType}
    @test A.alphabet == aa_alphabet

    A = Alignment(codonseqs; as_codons=true)
    @test A isa Alignment{PottsEvolver.Codon,PottsEvolver.IntType}
    @test A.alphabet == codon_alphabet

    A = Alignment(codonseqs; as_codons=false)
    @test A isa Alignment{Char,PottsEvolver.IntType}
    @test A.alphabet == aa_alphabet
end
