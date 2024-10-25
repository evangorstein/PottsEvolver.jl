module PottsEvolver

using ArgCheck
using BioSequenceMappings
using ProgressMeter
using Random
using StatsBase
using TreeTools
using UnPack

export read_fasta, symbols # from BioSequenceMappings
export read_tree # from TreeTools

import Base: ==, hash
import Base: copy, show, write
import Base: getindex, setindex!
import Base: iterate, length, eltype, size

import BioSequenceMappings: Alignment, to_string, hamming

# Default types for numerical quantities
const IntType = UInt8
const FloatType = Float64

include("codons.jl")
export codon_alphabet, aa_alphabet, nt_alphabet
export bases, genetic_code

include("sequences.jl")
export AbstractSequence, AASequence, CodonSequence
export Alignment
public translate

include("pottsgraph.jl")
export PottsGraph
export energy
public set_gauge!

include("sampling.jl")
export mcmc_sample, SamplingParameters
public get_init_sequence

include("sample_tree.jl")
public mcmc_sample_tree, pernode_alignment

include("IO.jl")
export read_graph, read_potts_graph

#=
- codons.jl: alphabets and genetic code
- sequences.jl: contain only a vector of Int (or two for CodonSequence).
  Conversion is done through alphabets
- IO.jl: for reading Potts models -- alignment/sequence IO is done through BioSequenceMappings
=#

end
