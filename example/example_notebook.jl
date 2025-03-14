### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a9825e4c-92e6-11ef-168d-0f2ba863e26d
begin
    using Pkg
    Pkg.activate(".")
	using Random
    using UnPack
    using PlutoUI
    using PottsEvolver
end

# ╔═╡ 992f54a2-c37a-47ce-9395-88dda801ecbc
begin
    using TreeTools # at this point it's just useful to have this
    using BioSequenceMappings
end

# ╔═╡ ab2d6cdb-b5be-438c-98c0-27bcb81f2f95
PlutoUI.TableOfContents()

# ╔═╡ 1c8c3060-3b4e-4407-a041-fc4f1c024332
md"# `PottsGraph`"

# ╔═╡ aed8bb56-d5e3-40ae-8133-11225486b62c
md"""
Potts models are contained in the `PottsGraph` type. They are relatively straightforward to manipulate
"""

# ╔═╡ ccfa1d4e-0bc7-4b10-9315-66d404f57574
# Read from file
potts = read_graph("parameters_PF00076.dat")

# ╔═╡ aa766af7-c050-4928-8c90-b1af5cbd3527
md"""
They essentially have three fields
- `h` and `J` contain fields and couplings, respectively of dimensions $(q \times L)$ and $(q \times q \times L \times L)$
- `alphabet` contains the mapping from integers (*e.g.* used to index `h` and `J`) and biological symbols (amino acids, codons...). It is an object of type `Alphabet`, coming from the `BioSequenceMappings` package. 
"""

# ╔═╡ 00979b31-56ae-45bc-891b-0c9d39fa0827
md"""
**Note**: In other places of the code, I kind of assume that the amino acid alphabet will be "$(prod(symbols(potts.alphabet)))". 

If your potts model does not use this, you can of course change `potts.alphabet` to be whatever, for instance `potts.alphabet = Alphabet("AC..."). However, I am not completely sure that this will go smoothly, please report any issue.
"""

# ╔═╡ 907c9922-b45c-4973-8bf3-7dcc8219a61f
[
    size(potts.h), # q*L
    size(potts.J), # q*q*L*L
    length(potts.alphabet), # q
    size(potts), # named tuple (; L, q)
]

# ╔═╡ 95d11aca-45ef-4720-8bc8-0812319c6854
(; L, q) = size(potts)

# ╔═╡ 58197d9e-3b67-40db-9852-6451cf4b254e
md"""
One can use the alphabet as a function to map from symbols to integers or the reverse. 
"""

# ╔═╡ 5045ef14-8547-43a2-a105-d4ec0962b04f
potts.alphabet

# ╔═╡ b8431fc3-d619-4ce3-bbf8-2e33e4963288
potts.alphabet(4) # convert integers to chars

# ╔═╡ 25e2f982-9a3c-4499-b2b9-dee74dc36951
potts.alphabet('H') # or the reverse

# ╔═╡ ed8559ec-2c67-4338-bcbb-20803700c24c
md"# Sampling an MCMC chain"

# ╔═╡ b1408a8a-f62d-46f2-ad05-a611324793ad
md"""
The main sampling function is `mcmc_sample`. See `?mcmc_sample` for details.
The following will sample a chain of length `M=10`, starting from a random codon sequence, and taking `Teq=100` steps between each sample. Finally, the output is translated back to amino acids
"""

# ╔═╡ 6576b8ba-762c-4e74-90b5-76f3725557d7
result = let
    M = 10
    parameters = SamplingParameters(; Teq=100, burnin=0)
    mcmc_sample(potts, M, parameters; init=:random_codon, translate_output=true)
end;

# ╔═╡ b7d03a71-eabc-43d9-97b4-8e388878718a
result

# ╔═╡ 7b174b8a-13bd-4534-bc0e-595149f9411f
md"""
The output contains the sample, the "time" of each sequence within the chain (*i.e.* how many steps away from the initial sequence), and other information.
"""

# ╔═╡ 7b455214-13e5-4646-a9c5-11e538475bca
begin
    @unpack sequences, tvals = result
    println("Time values: $tvals")
    sequences
end

# ╔═╡ f871ce6f-c80a-433e-a954-d34002d7b5c1
md"""
The sequences are in an `Alignment` object (from `BioSequenceMappings`). This essentially wraps the integer vectors in a `data` matrix, and the alphabet to convert them to biological symbols.

If the output was translated from codons to amino acids (as is the case here), the alphabet is automatically the same as the one from the input Potts model. Otherwise, it is `codon_alphabet`. 

We can use it to iterate over sequences, write the result to a fasta file, etc...
"""

# ╔═╡ 63fd2b58-4a71-4258-9eb4-96fb491bdd37
sequences[1] # the first sample (here also the initial sequence)

# ╔═╡ af3a221c-4998-43b0-9e2a-f1fb9da02536
[potts.alphabet(i) for i in sequences[1]]

# ╔═╡ 42e624ed-b38d-4908-a4fa-7bc635c78846
size(sequences) # length L x 10 sequences

# ╔═╡ 97aa28b3-012b-4614-9e2a-7706d9d48402
md"""
It's also easy to write / read alignments. 
"""

# ╔═╡ 0017a217-7cd6-4974-b936-30f93631396e
let
    write("example_alignment.fasta", sequences)
    X = read_fasta("example_alignment.fasta")
    all(x -> x[1] == x[2], zip(X, sequences))
end

# ╔═╡ 558f08ed-dfa0-4b18-8dc1-5c7c816b3534
md"""
### Parameters of the sampling
"""

# ╔═╡ 4a07a258-7a6f-4307-a759-ee9115020fea
md"""
They are provided by the `SamplingParameters` struct. 
The fields have the following meaning
- `Teq`: number of steps between samples (flips, not sweeps)
- `burnin`: number of steps between the initial sequence and the first sample
- `fraction_gap_step`: if using codon based sampling, fraction of metropolis gap steps. 
- `step_type`: only `:gibbs` for now
- `step_meaning`: what counts as a step when measuring "time" (`Teq`) between samples:
  + `:accepted` means only accepted steps count 
  + `:changed`: only steps leading to an actual state change count
  + `:proposed`: all steps count

`?SamplingParameters` may show more info.
"""

# ╔═╡ 8c2d4679-5c9a-46da-b053-2885fd4cf4a5
md"### Choosing the initial sequence"

# ╔═╡ 04fd2366-b9b5-4c45-ab0a-aee0d821fe1e
md"""
The initial sequence is provided to `mcmc_sample` either
- as a third positional argument
- as a kwarg `init`

In the first case, something of the type `AbstractSequence` should be used (see below). In the second case, the sequence is generated by calling `PottsEvolver.get_init_sequence(init, potts)`
the docstring of `PottsEvolver.get_init_sequence` will be helpful.

**Important**: the *type* of the initial sequence will determine the algorithm used for sampling:
- if it is a `CodonSequence`, then the MCMC is based on codons
- otherwise, it is just a standard MCMC
"""

# ╔═╡ aabe38f7-8f6e-4b66-949f-5da952193237
let
	Random.seed!(400)
    not_stop = filter(x -> x != 13 && x != 5, 1:21) # 
    s_not_valid = rand(not_stop, L) # Vector{Int} is not a valid input
    # Valid initial sequences
    aa_seq = AASequence(s_not_valid) # assumes the sequence is amino acids
    codon_seq = CodonSequence(s_not_valid) # assumes the sequence is codons
    num_seq = PottsEvolver.NumSequence(s_not_valid, 21) # assumes the sequence has no biological meaning

    # Using this for sampling
	parameters = SamplingParameters(; Teq=100, burnin=0)
	mcmc_sample(potts, 10, parameters; init=:random_aa) # operates on aa, but starting from a random sequence
    mcmc_sample(potts, 10, aa_seq, parameters) # use "normal" amino acid based mcmc 
    
    mcmc_sample(potts, 10, parameters; init=aa_seq) # this also works
    mcmc_sample(potts, 10, parameters; init=s_not_valid) # should work, but it's less clear
	res = mcmc_sample(potts, 10, codon_seq, parameters, translate_output=false) # use codon based mcmc
	res.sequences.alphabet
end

# ╔═╡ 65c82b80-cdfb-47f4-b398-989ec42ebbb3
let
    # Creating a random sequence is also easy
    rand_codon_seq = CodonSequence(L)
    # etc...
end

# ╔═╡ 58c13678-8cc7-45f2-9720-6b9b47596ccd
md"### Controlling the output"

# ╔═╡ 0fd48be5-b07f-4376-a915-c224990e2c97
md"""
The following arguments control the output of `run_mcmc`: 
- `alignment_output` (default: `true`): return an `Alignment` containing the sequences. If `false`, will return a vector of sequences. 
  The advantage of the alignment format is mainly that it is easy to convert to fasta, and that some convenient functions of `BioSequenceMappings` operate on alignment. 
- `translate_output` (default: `true`): if a `CodonSequence` was provided as input, the output would normally also consist of codon sequences (or a codon based alignment).  This flag causes the output to be translated back to amino acids. 
"""

# ╔═╡ 578b2fea-f418-4f3c-8423-5be10412c884
md"""
If you do not translate the output, here is what happens: 
"""

# ╔═╡ d9e5a4a1-7fab-4140-94dd-a4943bb16a0e
sequences_codon = let
    M = 10
    parameters = SamplingParameters(; Teq=100, burnin=0)
    mcmc_sample(potts, M, parameters; init=:random_codon).sequences
end;

# ╔═╡ 2ddc7737-8ae2-4725-91dd-5076c3bbd350
sequences_codon[1] # integers larger than 21 --> codons

# ╔═╡ 9bb03d89-8156-468d-a8e1-b7c8f753f9af
# Here are three ways to map this to amino acids
let
    alphabet = sequences_codon.alphabet # same as potts.alphabet
    s0 = sequences_codon[1]
    @info "Integer codons" s0 # Vector{Integer} standing for codons
    @info "Codons" alphabet(s0) # vector of PottsEvolver.Codon
    @info "Integer AA" map(x -> genetic_code(x), s0) # Vector{Integer} standing for amino acids
    @info "Char AA" map(x -> genetic_code(alphabet(x)), s0) # vector of Char, for amino acids
end

# ╔═╡ 434db5aa-92d4-45dc-abda-30757305aa5e
md"# Sampling on a tree"

# ╔═╡ 25215879-eb19-46a6-bbdf-15bd96e18000
md"""
*Note*: The package uses `TreeTools.jl` to manipulate trees, please have a look at its documentation if you want to handle trees. The `read_tree` function is re-exported from `TreeTools`. 

First, we read a small tree. It is a balanced binary tree with `8` leaves, and all branches of length `50`. 
"""

# ╔═╡ 4a522c90-0f74-4b35-a2c8-f64726f4875d
tree = read_tree("small_tree.nwk")

# ╔═╡ e1cfc82e-266c-48ac-85f3-c75f8c0cf7ca
for c in children(root(tree))
	println(branch_length(c))
end

# ╔═╡ c463db7d-cb3c-407b-b8ec-03ce04ade53c
md"""
Sampling on this tree is easy: just call `mcmc_sample` with the tree object instead of the number of sequences `M`. The number of steps on each branch will be equal to the length of the branch. 

!!! warn
	For now, branch lengths need to be integers (at least approximately). One thing I want to do is to add a parameter `:branch_length_interpretation` to change this. 
"""

# ╔═╡ b1aea1ba-c13d-4205-8761-3c6bf557da05
result_tree = let
    # Teq is useless here, but SamplingParameters complains if you don't set it
    parameters = SamplingParameters(; Teq=0, burnin=0)
    mcmc_sample(potts, tree, parameters; init=:random_codon)
end

# ╔═╡ 91103fb3-896b-420b-ba28-51871f725a11
md"""
The output contains
- the tree,
- an alignment of leaf sequences,
- an alignment of internal nodes sequences.

For each of these, the sequences are labeled using the label of tree nodes. 
"""

# ╔═╡ 77ed2514-772c-45be-8b04-f4593d80a462
let
    # Find a sequence given a leaf label
    interesting_leaf = label(collect(leaves(tree))[5])
    @info "Let's look at leaf $(interesting_leaf)"
    seq = find_sequence(interesting_leaf, result_tree.leaf_sequences)[2]
    @info "Corresponding sequence: $(result_tree.leaf_sequences.alphabet(seq))"

	interesting_node = label(collect(internals(tree))[6])
    @info "Let's look at node $(interesting_node)"
    seq = find_sequence(interesting_node, result_tree.internal_sequences)[2]
    @info "Corresponding sequence: $(result_tree.internal_sequences.alphabet(seq))"
	
end

# ╔═╡ e7e26419-29d1-4444-a248-5979c4ba755f
md"# Practical example"

# ╔═╡ eeab33e7-a477-4cac-824b-4b084fccb7fa
md"""
Using a fixed tree, let's simulate evolution many times. 
We want to start from equilibrated sequences and sample `K` times. 
"""

# ╔═╡ a0330814-385f-4296-aa70-7d111267bb1f
let
    K = 10
    # Long mcmc chain to sample K sequences
    parameters = SamplingParameters(; step_meaning=:accepted, Teq=100 * L, burnin=1000 * L)
    # no translate because we want to sample on tree using the codon algorithm
    # no alignment output because we want the CodonSequence at the root
    init_sequences =
        mcmc_sample(
            potts,
            K,
            parameters;
            init=:random_codon,
            translate_output=false,
            alignment_output=false,
        ).sequences
    @info typeof(init_sequences)

    # For each initial sequence, sample on the tree
    # save result to some directory
    savedir = mkpath("_scrap")
    parameters = SamplingParameters(; step_meaning=:accepted, Teq=0)
    for (i, init) in enumerate(init_sequences)
        results = mcmc_sample(potts, tree, parameters; init)
        @unpack leaf_sequences, internal_sequences = results
        write(joinpath(savedir, "leaf_sequences_$i.fasta"), leaf_sequences)
        write(joinpath(savedir, "internal_sequences_$i.fasta"), internal_sequences)
    end
    # write the tree for reference
    write(joinpath(savedir, "tree.nwk"), tree)
end

# ╔═╡ Cell order:
# ╠═a9825e4c-92e6-11ef-168d-0f2ba863e26d
# ╠═ab2d6cdb-b5be-438c-98c0-27bcb81f2f95
# ╟─1c8c3060-3b4e-4407-a041-fc4f1c024332
# ╟─aed8bb56-d5e3-40ae-8133-11225486b62c
# ╠═ccfa1d4e-0bc7-4b10-9315-66d404f57574
# ╟─aa766af7-c050-4928-8c90-b1af5cbd3527
# ╟─00979b31-56ae-45bc-891b-0c9d39fa0827
# ╠═907c9922-b45c-4973-8bf3-7dcc8219a61f
# ╠═95d11aca-45ef-4720-8bc8-0812319c6854
# ╟─58197d9e-3b67-40db-9852-6451cf4b254e
# ╠═5045ef14-8547-43a2-a105-d4ec0962b04f
# ╠═b8431fc3-d619-4ce3-bbf8-2e33e4963288
# ╠═25e2f982-9a3c-4499-b2b9-dee74dc36951
# ╟─ed8559ec-2c67-4338-bcbb-20803700c24c
# ╟─b1408a8a-f62d-46f2-ad05-a611324793ad
# ╠═6576b8ba-762c-4e74-90b5-76f3725557d7
# ╠═b7d03a71-eabc-43d9-97b4-8e388878718a
# ╟─7b174b8a-13bd-4534-bc0e-595149f9411f
# ╠═7b455214-13e5-4646-a9c5-11e538475bca
# ╟─f871ce6f-c80a-433e-a954-d34002d7b5c1
# ╠═63fd2b58-4a71-4258-9eb4-96fb491bdd37
# ╠═af3a221c-4998-43b0-9e2a-f1fb9da02536
# ╠═42e624ed-b38d-4908-a4fa-7bc635c78846
# ╟─97aa28b3-012b-4614-9e2a-7706d9d48402
# ╠═0017a217-7cd6-4974-b936-30f93631396e
# ╟─558f08ed-dfa0-4b18-8dc1-5c7c816b3534
# ╟─4a07a258-7a6f-4307-a759-ee9115020fea
# ╟─8c2d4679-5c9a-46da-b053-2885fd4cf4a5
# ╟─04fd2366-b9b5-4c45-ab0a-aee0d821fe1e
# ╠═aabe38f7-8f6e-4b66-949f-5da952193237
# ╠═65c82b80-cdfb-47f4-b398-989ec42ebbb3
# ╟─58c13678-8cc7-45f2-9720-6b9b47596ccd
# ╟─0fd48be5-b07f-4376-a915-c224990e2c97
# ╟─578b2fea-f418-4f3c-8423-5be10412c884
# ╠═d9e5a4a1-7fab-4140-94dd-a4943bb16a0e
# ╠═2ddc7737-8ae2-4725-91dd-5076c3bbd350
# ╠═9bb03d89-8156-468d-a8e1-b7c8f753f9af
# ╟─434db5aa-92d4-45dc-abda-30757305aa5e
# ╟─25215879-eb19-46a6-bbdf-15bd96e18000
# ╠═4a522c90-0f74-4b35-a2c8-f64726f4875d
# ╠═e1cfc82e-266c-48ac-85f3-c75f8c0cf7ca
# ╟─c463db7d-cb3c-407b-b8ec-03ce04ade53c
# ╠═b1aea1ba-c13d-4205-8761-3c6bf557da05
# ╟─91103fb3-896b-420b-ba28-51871f725a11
# ╠═992f54a2-c37a-47ce-9395-88dda801ecbc
# ╠═77ed2514-772c-45be-8b04-f4593d80a462
# ╟─e7e26419-29d1-4444-a248-5979c4ba755f
# ╟─eeab33e7-a477-4cac-824b-4b084fccb7fa
# ╠═a0330814-385f-4296-aa70-7d111267bb1f
