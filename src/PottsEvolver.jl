module PottsEvolver

using BioSequenceMappings
using ProgressMeter
using Random
using StatsBase
using TreeTools

export hamming, symbols # from BioSequenceMappings

import Base: ==, hash
import Base: copy, show, size, write
import Base: getindex, setindex!

import BioSequenceMappings: Alignment, to_string, hamming

include("pottsgraph.jl")
export PottsGraph

include("codons.jl")
export codon_alphabet, aa_alphabet, nt_alphabet
export bases, genetic_code
export isstop, isgap

include("sequences.jl")
export AbstractSequence, AASequence, CodonSequence
export Alignment
export translate

include("sampling.jl")
export mcmc_sample, SamplingParameters

include("IO.jl")
export read_graph, read_potts_graph

#=
- codons.jl: alphabets and genetic code
- sequences.jl: contain only a vector of Int (or two for CodonSequence).
  Conversion is done through alphabets


=#

end
