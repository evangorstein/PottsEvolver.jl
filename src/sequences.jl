abstract type AbstractSequence end

Base.length(s::AbstractSequence) = length(s.seq)
Base.getindex(s::AbstractSequence, i) = getindex(s.seq, i)
Base.setindex!(s::AbstractSequence, x, i) = setindex!(s.seq, x, i)
function Base.:(==)(x::T, y::T) where T <: AbstractSequence
    return all(p -> getproperty(x, p) == getproperty(y, p), propertynames(x))
end
function Base.hash(x::AbstractSequence, h::UInt)
    return hash(x.seq, h)
end

# the two functions below: used in Alignment
sequence(x::AbstractSequence; kwargs...) = x.seq
_sequence_alphabet(::Type{<:AbstractSequence}; kwargs...) = nothing

#=
Other methods that a subtype should implement
- copy
=#

#====================================#
############# AASequence #############
#====================================#

mutable struct AASequence{T<:Integer} <: AbstractSequence
    seq::Vector{T}
    function AASequence(x::Vector{T}) where T
        q = length(aa_alphabet)
        @assert all(<=(q), x) "AA are represented by `(1..$(q))` integers. Instead, $x"
        return new{T}(x)
    end
end

Base.copy(s::AASequence) = AASequence(copy(s.seq))
AASequence(L::Integer; T = IntType) = AASequence(rand(T(1):T(length(aa_alphabet)), L))

_sequence_alphabet(::Type{<:AASequence}; kwargs...) = aa_alphabet

#=============================================#
################ CodonSequence ################
#=============================================#

mutable struct CodonSequence{T<:Integer} <: AbstractSequence
    seq::Vector{T} # the codons
    aaseq::Vector{T} # the translation
    function CodonSequence(seq::Vector{T}, aaseq::Vector{T}) where T
        qc = length(codon_alphabet)
        qaa = length(aa_alphabet)
        @assert all(<=(qc), seq) """
            Codons are represented by `(1..$(qc))` integers. Instead $seq
        """
        @assert all(<=(qaa), aaseq) """
            AA are represented by `(1..$(qaa))` integers. Instead, $aaseq
        """
        any(isstop, seq) && @warn "Sequence contains stop codon"
        @assert all(x -> genetic_code(x[1]) == x[2], zip(seq, aaseq)) """
            Codon and amino acid sequences do not match. Got $seq and $aaseq
        """
        return new{T}(seq, aaseq)
    end
end

"""
    CodonSequence(seq::Vector{Integer}; source=:aa)

Build a `CodonSequence` from `seq`:
- if `source==:codon`, `seq` is interpreted as representing codons (see `codon_alphabet`);
- if `source==:aa`, `seq` is interpreted as representing amino acids (see `aa_alphabet`);
  matching codons are randomly chosen using the `PottsEvolver.reverse_code_rand` method.
"""
function CodonSequence(seq::Vector{T}; source=:codon) where T <: Integer
    return if source == :aa
        CodonSequence(map(reverse_code_rand, seq), convert(Vector{T}, seq))
    elseif source == :codon
        aaseq = map(genetic_code, seq)
        any(isnothing, aaseq) && error("""
            Cannot build `CodonSequence` from input that contains stop codon.
            Input sequence was $seq."""
        )
        CodonSequence(convert(Vector{T}, seq), aaseq)
    else
        error("Unknown `source` $(source). Use `:aa` or `:codon`.")
    end
end
"""
    CodonSequence(L::Int; source=:aa, T)

Sample `L` states at random of the type of `source` (`:aa` or `:codon`).
Underlying integer type is `T`.
"""
function CodonSequence(L::Int; source=:aa, T = IntType)
    return CodonSequence(rand(T(1):T(length(aa_alphabet)), L); source)
end

function Base.setindex!(s::CodonSequence, x::Integer, i)
    isstop(x) && @warn "Introducing stop codon in sequence"
    setindex!(s.aaseq, genetic_code(x), i)
    return setindex!(s.seq, x, i)
end
Base.copy(s::CodonSequence) = CodonSequence(copy(s.seq), copy(s.aaseq))

translate(s::CodonSequence) = AASequence(s.aaseq)
function sequence(x::CodonSequence; as_aa=false, as_codons=true)
    if !xor(as_aa, as_codons)
        error("Expected either `as_aa` or `as_codons`. Instead $as_aa $as_codons")
    end
    as_aa ? x.aaseq : x.seq
end

function _sequence_alphabet(::Type{<:CodonSequence}; as_aa=false, as_codons=true)
    if !xor(as_aa, as_codons)
        error("Expected either `as_aa` or `as_codons`. Instead $as_aa $as_codons")
    end
    return as_aa ? aa_alphabet : codon_alphabet
end

#============================================================#
##################### Numerical sequence #####################
#============================================================#

@kwdef mutable struct NumSequence{T<:Integer} <: AbstractSequence
    seq::Vector{T}
    q::T = maximum(seq)
    function NumSequence(seq::AbstractVector{T}, q) where T
        @assert all(<=(q), seq) "All entries in `seq` must be smaller than $q. Instead $seq"
        return new{T}(seq, q)
    end
end
NumSequence(L::Integer, q::Integer; T = IntType) = NumSequence(rand(T(1):T(q), L), T(q))



Base.copy(x::NumSequence) = NumSequence(copy(x.seq), x.q)

"""
    NumSequence(L, q; T)

Construct a random sequence of integers of length `L` using integers `1:q`.
The integer type can be set using `T`.
"""

#===========================================================================#
########################## Converting to Alignment ##########################
#===========================================================================#

"""
    Alignment(sequences; alphabet, names, as_aa, as_codons)

Construct a `BioSequenceMappings.Alignment` from a set of sequences.
If the sequences are `AASequence`, `alphabet` defaults to `aa_alphabet`.
If they are `CodonSequence`, `as_aa` or `as_codon` can be used to decide whether the
alignment should store codons or amino acids. `alphabet` can be determined automatically
from this.
"""
function Alignment(
    S::AbstractVector{T};
    names = nothing, as_aa=true, as_codons=false,
    alphabet = _sequence_alphabet(eltype(S); as_aa, as_codons),
) where T <: AbstractSequence

    # Checks
    if !xor(as_aa, as_codons)
        error("Expected either `as_aa` or `as_codons`. Instead $as_aa $as_codons")
    end
    if !allequal(Iterators.map(length, S))
        error("Sequences do not have the same length")
    end
    if !isnothing(names) && length(names) != length(S)
        error("Got $(length(names)) but $(length(S)) sequences.")
    end

    data = hcat([sequence(s; as_aa, as_codons) for s in S]...)

    return Alignment(data, alphabet; names)
end





# function Alignment(
#     S::AbstractVector{CodonSequence};
#     as_aa=true, as_codons=false, alphabet = nothing, names = nothing,
# )
#     # Checks
#     if !xor(as_aa, as_codons)
#         error("Expected either `as_aa` or `as_codons`. Instead $as_aa $as_codons")
#     end
#     if !allequal(Iterators.map(length, S))
#         error("Sequences do not have the same length")
#     end
#     if !isnothing(names) && length(names) != length(S)
#         error("Got $(length(names)) but $(length(S)) sequences.")
#     end

#     if isnothing(alphabet)
#         alphabet = as_aa ? aa_alphabet : codon_alphabet
#     end

#     data = Matrix{IntType}(undef, length(first(S)), length(S))
#     for (m, s) in enumerate(S)
#         data[:, m] .= as_aa ? s.aaseq : s.seq
#     end

#     return Alignment(data, alphabet; names)
# end


#==================#
####### Misc #######
#==================#

function BioSequenceMappings.hamming(x::AbstractSequence, y::AbstractSequence; kwargs...)
    return hamming(x.seq, y.seq; kwargs...)
end
function BioSequenceMappings.hamming(
    x::CodonSequence, y::CodonSequence;
    source = :codon, kwargs...,
)
    return if source == :codon
        hamming(x.seq, y.seq; kwargs...)
    elseif source == :aa
        hamming(x.aaseq, y.aaseq; kwargs...)
    else
        error("Valid `source` values: `:codon` or `:aa`. Instead $source")
    end
end
