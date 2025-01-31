#============================================================#
##################### SamplingParameters #####################
#============================================================#

const VALID_STEP_TYPES = [:gibbs, :metropolis]
const VALID_STEP_MEANINGS = [:proposed, :accepted, :changed]

"""
```
BranchLengthMeaning

    type::Symbol
    length::Symbol
```

Define how the branch lengths should be interpreted when sampling on a tree.
This is done in two consecutive steps.
1. Convert the branch lengths to mcmc steps.
    If `type==:sweep`, multiply the branch length by the length of the sequence.
    Else, if `type==:step`, leave it as is.
2. Convert the branch length to an integer value which will be the number of mcmc steps.
  - if `length==:exact`, the branch length is expected to be an integer (`rtol=1e-6`).
  - if `length==:round`, round it to the closest `Int`.
  - if `length==:poisson`, sample a poisson distribution with the branch length as mean.
"""
@kwdef struct BranchLengthMeaning
    type::Symbol
    length::Symbol
    function BranchLengthMeaning(type, length)
        valid_types = (:step, :sweep)
        valid_lengths = (:exact, :round, :poisson)
        @argcheck type in valid_types """
        `type` field should be in `$(valid_types)`. Instead `:$type`.
        """
        @argcheck length in valid_lengths """
        `length` field should be in `$(valid_lengths)`. Instead `:$length`.
        """
        return new(type, length)
    end
end

function steps_from_branchlength(τ::Real, p::BranchLengthMeaning, L::Int)
    if p.type == :sweep
        τ *= L
    end

    rtol = 1e-6
    if p.length == :exact
        τi = round(Int, τ)
        if !isapprox(τ, τi; rtol)
            msg = """
            Expected tree with approx integer branch lengths (`rtol=$rtol`). Instead $τ.
            """
            throw(ArgumentError(msg))
        end
        return τi
    elseif p.length == :round
        return round(Int, τ)
    elseif p.length == :poisson
        return pois_rand(τ)
    end
    @assert false "This state should never be reached"
end

"""
    mutable struct SamplingParameters

Construct using keyword arguments:
```
step_type::Symbol = :gibbs
step_meaning::Symbol = :accepted
Teq::Int
burnin::Int = 5*Teq
fraction_gap_step::Float64 = 0.9
branchlength_meaning::BranchLengthMeaning
```

- `Teq` is measured in swaps: attempted (or accepted) change of one sequence position.
- `burnin`: number of steps starting from the initial sequence before the first `Teq` are
    made.
- `step_meaning` can take three values
    - `:proposed`: all mcmc steps count towards equilibration
    - `:accepted`: only accepted steps count (all steps for non-codon Gibbs)
    - `:changed`: only steps that lead to a change count (Gibbs can resample the same state)
    *Note*: Gibbs steps for codons are more complicated, since they involve the possibility
    Metropolis step for gaps, which can be rejected.
- `branchlength_meaning` is only useful if you sample along branches of a tree.
  See `?BranchLengthMeaning` for information.
"""
@kwdef mutable struct SamplingParameters
    step_type::Symbol = :gibbs
    step_meaning::Symbol = :accepted
    Teq::Int
    burnin::Int = 5 * Teq
    fraction_gap_step::Float64 = 0.9
    branchlength_meaning::BranchLengthMeaning = BranchLengthMeaning(:step, :exact)
    function SamplingParameters(
        step_type, step_meaning, Teq, burnin, fraction_gap_step, branchlength_meaning
    )
        step_meaning = try
            Symbol(step_meaning)
        catch err
            @error "Invalid `step_meaning` $step_meaning"
            throw(err)
        end
        @argcheck step_meaning in VALID_STEP_MEANINGS """
                `step_meaning` should be in $VALID_STEP_MEANINGS.
                Instead $(step_meaning).
            """
        @argcheck step_type in VALID_STEP_TYPES """
                `step_type` should be in $VALID_STEP_TYPES.
                Instead $(step_type).
            """
        return new(
            step_type, step_meaning, Teq, burnin, fraction_gap_step, branchlength_meaning
        )
    end
end

#=======================================#
############## mcmc_steps! ##############
#=======================================#

"""
    mcmc_steps!(
        s::AbstractSequence, g, num_steps, p::SamplingParameters;
        gibbs_holder, kwargs...
    )

Perform `num_steps` MCMC steps starting from sequence `s` and using graph `g`.
The step type (`:gibbs`, `:metropolis`) and the interpretation of `num_steps`
    (`:changed`, `:accepted`, `:proposed`) is set using `p` (see `?SamplingParameters`).
Modifies the input sequence `s` and returns it.
"""
function mcmc_steps!(
    sequence::AbstractSequence,
    g::PottsGraph,
    num_steps::Integer,
    p::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    gibbs_holder=get_gibbs_holder(sequence),
    verbose=false,
)
    @argcheck length(sequence) == size(g).L """
    Size of sequence and model do not match: $(length(sequence)) and $(size(g).L)
    """

    step_func! = if p.step_type == :gibbs
        gibbs_step!
    elseif p.step_type == :metropolis
        metropolis_step!
    else
        error("Unknown `step_type` $(p.step_type)")
    end

    proposed = 0
    accepted = 0
    if p.step_meaning == :proposed
        for _ in 1:num_steps
            step_func!(sequence, g, p, gibbs_holder; rng)
            proposed += 1
            accepted += 1
        end
    else
        min_acceptance_rate = 0.001
        max_tries = num_steps / min_acceptance_rate
        while accepted < num_steps && proposed < max_tries
            step_result = step_func!(sequence, g, p, gibbs_holder; rng)
            if p.step_meaning == :accepted && step_result.accepted
                accepted += 1
            elseif p.step_meaning == :changed && step_result.changed
                accepted += 1
            end
            proposed += 1
        end
        proposed >= max_tries && max_tries > 0 && @warn """
            $max_tries steps attempted with only $accepted accepted. Giving up.
            """
    end

    return sequence, proposed, accepted
end

#=
MCMC step functions
the return value should be a named tuple of the form
`(i, new_state, accepted, changed)`, plus extra info if needed
`accepted` says whether the step was accepted (always true for non-codon gibbs)
`change` is specific to codon sequences, and says whether the amino acid state was changed.

the input arguments should be
`(sequence, graph, params, holder)` where holder is only really used for gibbs
=#

#===============================================================#
###################### CodonSequence steps ######################
#===============================================================#

function gibbs_step!(
    s::CodonSequence, g::PottsGraph, params::SamplingParameters, gibbs_holder; kwargs...
)
    p = params.fraction_gap_step
    return if rand() < p
        gap_metropolis_step!(s, g; kwargs...)
    else
        aa_gibbs_step!(s, g, gibbs_holder; kwargs...)
    end
end

function metropolis_step!(s::CodonSequence, g::PottsGraph, p::SamplingParameters; kwargs...)
    return error("Not implemented")
end

"""
    aa_gibbs_step!(s::CodonSequence, g::PottsGraph; kwargs...)

Change one coding codon in `s` into another coding codon.
"""
function aa_gibbs_step!(
    s::CodonSequence, g::PottsGraph, gibbs_holder; rng=Random.GLOBAL_RNG
)
    p = gibbs_holder
    # new_codons and new_aas are arrays but should _NOT_ be mutated
    i, b, new_codons, new_aas = pick_aa_mutation(s; rng)
    aaref = s.aaseq[i]

    # Constructing the gibbs field p
    for (a, aanew) in enumerate(new_aas)
        if aanew == aaref
            p[a] = 0.0 # energy for now
        else
            # ΔE is minus the energy --> softmax(ΔE) is the Boltzmann distribution
            # ~ high ΔE is more probable
            ΔE = -aa_degeneracy(aanew) + aa_degeneracy(aaref) # log-degeneracy of gen code
            ΔE += g.h[aanew, i] - g.h[aaref, i]
            for j in 1:length(s)
                if j != i
                    ΔE += g.J[aanew, s.aaseq[j], i, j] - g.J[aaref, s.aaseq[j], i, j]
                end
            end
            p[a] = g.β * ΔE
        end
    end

    # Sampling from p
    if length(new_aas) < length(p)
        p[(length(new_aas) + 1):end] .= -Inf
    end
    softmax!(p)
    c = sample_from_weights(p)
    s[i] = new_codons[c]
    s.aaseq[i] = new_aas[c]
    changed = (new_aas[c] != aaref)
    accepted = true

    return (i=i, new_state=s[i], accepted=accepted, changed=changed)
end
function aa_gibbs_step!(
    s::CodonSequence, g::PottsGraph{T}; gibbs_holder=zeros(T, 4), kwargs...
) where {T}
    return aa_gibbs_step!(s, g, gibbs_holder; kwargs...)
end

function gap_metropolis_step!(
    s::CodonSequence, g::PottsGraph; rng=Random.GLOBAL_RNG, kwargs...
)
    i = rand(1:length(s))
    codon_ref = s.seq[i]
    aa_ref = s.aaseq[i]
    # di Bari et. al.
    β = 1 / n_aa_codons # aa to gap transition
    if codon_ref == gap_codon_index
        # Pick any non gap codon and try to mutate to it
        codon_new = rand(coding_codons)
        aa_new = genetic_code(codon_new)
        ΔE = -aa_degeneracy(aa_new) + aa_degeneracy(aa_ref)
        ΔE += g.h[aa_new, i] - g.h[aa_ref, i]
        for j in 1:length(s)
            if j != i
                ΔE += g.J[aa_new, s.aaseq[j], i, j] - g.J[aa_ref, s.aaseq[j], i, j]
            end
        end

        return _metropolis_step!(s, ΔE, i, codon_new, aa_new, aa_ref)
    else
        # Try to replace s[i] by the gap codon
        if rand() < 1 - β
            # with probability 1-β do nothing
            # return value below should match _metropolis_step!(::CodonSequence, ...)
            (i=i, new_state=s[i], accepted=false, changed=false)
        else
            # Metropolis move: replace s[i] with gap codon
            codon_new = gap_codon_index # gap_codon_index defined in codons.jl
            aa_new = genetic_code(codon_new)
            ΔE = aa_degeneracy(aa_ref) # degeneracy of gap state (aa_new) is 0
            ΔE += g.h[aa_new, i] - g.h[aa_ref, i]
            for j in 1:length(s)
                if j != i
                    ΔE += g.J[aa_new, s.aaseq[j], i, j] - g.J[aa_ref, s.aaseq[j], i, j]
                end
            end

            return _metropolis_step!(s, ΔE, i, codon_new, aa_new, aa_ref)
        end
    end
end

function _metropolis_step!(s::CodonSequence, ΔE, i, codon_new, aa_new, aa_ref)
    # Check for acceptance based on ΔE
    # if accepted, change s.seq and s.aaseq using codon_new/aa_new
    # otherwise, do nothing
    return if ΔE >= 0 || rand() < exp(ΔE)
        s[i] = codon_new
        s.aaseq[i] = aa_new
        (i=i, new_state=s[i], accepted=true, changed=aa_new == aa_ref)
    else
        (i=i, new_state=s[i], accepted=false, changed=false)
    end
end

#===================================================#
################## Non-codon steps ##################
#===================================================#

function gibbs_step!(
    s::AbstractSequence,
    g::PottsGraph,
    p::SamplingParameters,
    gibbs_holder;
    rng=Random.GLOBAL_RNG,
)
    p = gibbs_holder
    q = length(gibbs_holder)

    i = rand(1:length(s))
    a_ref = s[i]
    for a in 1:q
        if a == a_ref
            p[a] = 0.0
        else
            ΔE = g.h[a, i] - g.h[a_ref, i]
            for j in 1:length(s)
                if j != i
                    ΔE += g.J[a, s[j], i, j] - g.J[a_ref, s[j], i, j]
                end
            end
            p[a] = g.β * ΔE
        end
    end

    softmax!(p)
    c = wsample(p)
    s[i] = wsample(p)
    changed = (a_ref != s[i])
    return (; i, new_state=s[i], accepted=true, changed)
end

#=====================#
######## Utils ########
#=====================#

"""
    pick_aa_mutation(s::CodonSequence)

Pick a valid codon to mutate in `s`.
Return the position `i`, the base `b`, new potential codons and aas.
Gap codons cannot be mutated using this.
"""
function pick_aa_mutation(s::CodonSequence; rng=Random.GLOBAL_RNG)
    max_cnt = 1000
    cnt = 1
    # Pick i and b, and see if we can mutate the codon (if it's a gap we can't)
    @label again
    cnt += 1
    i = rand(1:length(s))
    b = rand(IntType(1):IntType(3))
    codon = s[i]
    new_codons, new_aas = accessible_codons(codon, b)
    if cnt < max_cnt && (isnothing(new_codons) || isnothing(new_aas))
        @goto again
    else
        return i, b, new_codons, new_aas
    end
end

function softmax!(X)
    # thanks chatGPT!
    # Compute the maximum value in X to ensure numerical stability
    max_val = maximum(X)
    Z = 0.0
    @inbounds for i in eachindex(X)
        X[i] -= max_val
        X[i] = exp(X[i])
        Z += X[i]
    end
    @inbounds for i in eachindex(X)
        X[i] /= Z
    end
    return X
end

function sample_from_weights(W)
    # !!! Assumes W is normalized !!!
    x = rand()
    z = 0.0
    for (i, w) in enumerate(W)
        z += w
        if x < z
            return i
        end
    end
    return length(W)
end

get_gibbs_holder(::CodonSequence, T=FloatType) = zeros(T, 4)
get_gibbs_holder(::AASequence, T=FloatType) = zeros(T, 21)
get_gibbs_holder(s::NumSequence, T=FloatType) = zeros(T, s.q)
