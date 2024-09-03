#==================================================================#
####################### Alphabets and Codons #######################
#==================================================================#

const nucleotides = ['A', 'C', 'G', 'T']

const IntType = UInt8
"""
    aa_alphabet
    nt_alphabet
    codon_alphabet

Default alphabets for `PottsEvolver`.
"""
const aa_alphabet = Alphabet(Alphabet(:aa).characters, IntType; default_char=nothing)

"""
    aa_alphabet
    nt_alphabet
    codon_alphabet

Default alphabets for `PottsEvolver`.
"""
const nt_alphabet = Alphabet(nucleotides, IntType)

@kwdef struct Codon
    b1::Char
    b2::Char
    b3::Char
end
"""
    bases(codon)

Return iterator on the bases of `codon`.
"""
bases(codon::Codon) = Iterators.map(i -> getfield(codon, i), 1:3)

"""
    aa_alphabet
    nt_alphabet
    codon_alphabet

Default alphabets for `PottsEvolver`.
"""
const codon_alphabet = let
    nt = nucleotides
    C = map(x -> Codon(x...), Iterators.product(nt, nt, nt)) |> vec
    pushfirst!(C, Codon('-', '-', '-'))
    Alphabet(C, IntType)
end

#==========================================#
############### Genetic code ###############
#==========================================#

const aa_order = ['K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F']
const _genetic_code_struct = let
    code = Dict{Codon, Char}()
    i = 1
    for a in nucleotides, b in nucleotides, c in nucleotides
        code[Codon(a,b,c)] = aa_order[i]
        i += 1
    end
    code[Codon('-', '-', '-')] = '-'
    code
end
const _genetic_code_integers = let
    code = Dict{IntType, Union{Nothing, IntType}}()
    for (codon, aa) in _genetic_code_struct
        code[codon_alphabet(codon)] = aa == '*' ? nothing : aa_alphabet(aa)
    end
    code
end
"""
    genetic_code(x::Integer)

Translate the `i`th codon and return the index of the corresponding amino acid, using
the default `aa_alphabet`
"""
function genetic_code(codon::T) where T <: Integer
    aa = _genetic_code_integers[codon]
    isnothing(aa) ? aa : T(aa)
end
"""
    genetic_code(c::Codon)

Translate `c` and return the amino acid as a `Char`.
"""
genetic_code(codon::Codon) = _genetic_code_struct[codon]


"""
    reverse_code(aa)

Return the set of codons coding for `aa`.
"""
function reverse_code(aa::T) where T <: Integer
    return T.(findall(c -> genetic_code(c)==aa_alphabet(aa), codon_alphabet.index_to_char))
end
reverse_code(aa::AbstractChar) = map(codon_alphabet, reverse_code(aa_alphabet(aa)))
"""
    reverse_code_rand(aa)

Return a random codon coding for `aa`
"""
function reverse_code_rand(aa::Integer)
    codons = reverse_code(aa)
    # @info codons
    if isempty(codons)
        error("No codon corresponds to amino acid $aa ($aa_alphabet(aa))")
    end
    return rand(codons)
end
reverse_code_rand(aa::AbstractChar) = aa_alphabet(aa) |> reverse_code_rand |> codon_alphabet


#========================================================================#
######################### Codon helper functions #########################
#========================================================================#


function Base.show(io::IO, c::Codon)
    aa, iaa = if !isstop(c)
        aa = genetic_code(c)
        iaa = aa_alphabet(aa)
        aa, iaa
    else
        aa = '*'
        iaa = "STOP"
        aa, iaa
    end
    # print(io, "Codon $(codon_alphabet(c)): \"$(c.b1)$(c.b2)$(c.b3)\" --> $aa($iaa)")
    print(io, "\"$(c.b1)$(c.b2)$(c.b3)\"")
end
function Base.show(io::IO, x::MIME"text/plain", c::Codon)
    aa, iaa = if !isstop(c)
        aa = genetic_code(c)
        iaa = aa_alphabet(aa)
        aa, iaa
    else
        aa = '*'
        iaa = "STOP"
        aa, iaa
    end
    println(io, "Codon $(codon_alphabet(c)): \"$(c.b1)$(c.b2)$(c.b3)\" --> $aa($iaa)")
end

# Only all gaps or all nt codons are valid
# Anything else means a frameshift and I do not deal with that here
isgap(c::Codon) = all(==('-'), bases(c))
isgap(i::Integer) = isgap(codon_alphabet(i))
function isvalid(c::Codon)
    return if isgap(c)
        true
    elseif all(in(nucleotides), bases(c))
        true
    else
        false
    end
end
isstop(c::Codon) = genetic_code(c) == '*'
isstop(i::Integer) = isstop(codon_alphabet(i))

function to_string(s::AbstractVector{<:Integer}, alphabet::Alphabet{Codon, <:Integer})
    return prod(x -> prod(bases(alphabet(x))), s)
end

#===========================================================================#
########################## Codon accessibility map ##########################
#===========================================================================#

#=
For a given codon `c` and a given position `i` inside the codon,
what other codons are accessible by one mutation?

`_codon_access_map`: let `c::Int` be a codon and `i ∈ [1,3]` a position. `codon_access_map[c,i]` returns a tuple with:
    - the list of codons accessible by mutating `c` at position `i`, as integers. Stop codons or invalid codons are filtered out.
    - the list of corresponding amino acids, again as integers
!!! Only nucleotide mutations are considered. For this reason the gap codon never appears in this dictionary.
!!! I consider here that a codon is always accessible from itself! *i.e.* `c` will appear in the list
=#

const _codon_access_map = let
    M = Dict{Tuple{IntType,IntType}, Tuple{Vector{IntType}, Vector{IntType}}}()
    for c in 1:length(codon_alphabet), i in 1:3
        codon = codon_alphabet(c)
        if !isgap(codon) && isvalid(codon)
            accessible_codons = map(nucleotides) do a
                nts = [codon.b1, codon.b2, codon.b3]
                nts[i] = a
                codon_alphabet(Codon(nts...))
            end
            filter!(!isstop, accessible_codons)
            M[c,i] = (accessible_codons, map(genetic_code, accessible_codons))
        end
    end
    M
end
"""
    accessible_codons(codon, b::Integer)

Return all codons/amino-acids accessible by mutating `codon` at base `b`.
Value returned is a `Tuple` whose first/second elements represent codons/amino-acids.
"""
function accessible_codons(codon::T, b::Integer) where T<:Integer
    return get(_codon_access_map, (codon, b), (nothing, nothing))
end
function accessible_codons(codon::Codon, b::Integer)
    Cs = get(_codon_access_map, (codon_alphabet(codon), b), nothing)
    return if isnothing(Cs)
        Codon[], Char[]
    else
        map(codon_alphabet, Cs[1]), map(aa_alphabet, Cs[2])
    end
end
