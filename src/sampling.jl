const VALID_STEP_TYPES = [:gibbs, :metropolis]
const VALID_STEP_MEANINGS = [:proposed, :accepted, :changed]

"""
    mutable struct SamplingParameters

Construct using keyword arguments:
```
step_type::Symbol = :gibbs
step_meaning::Symbol = :accepted
Teq::Int
burnin::Int = 5*Teq
fraction_gap_step::Float64 = 0.9
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
"""
@kwdef mutable struct SamplingParameters
    step_type::Symbol = :gibbs
    step_meaning::Symbol = :accepted
    Teq::Int
    burnin::Int = 5 * Teq
    fraction_gap_step::Float64 = 0.9
    function SamplingParameters(step_type, step_meaning, Teq, burnin, fraction_gap_step)
        @argcheck step_meaning in VALID_STEP_MEANINGS """
                `step_meaning` should be in $VALID_STEP_MEANINGS.
                Instead $(step_meaning).
            """
        @argcheck step_type in VALID_STEP_TYPES """
                `step_type` should be in $VALID_STEP_TYPES.
                Instead $(step_type).
            """
        return new(step_type, step_meaning, Teq, burnin, fraction_gap_step)
    end
end

"""
    mcmc_sample(
        g::PottsGraph, M::Integer, s0::AbstractSequence, params::SamplingParameters;
        rng=Random.GLOBAL_RNG, verbose=0, progress_meter=true, alignment_output=true,
    )

Sample `g` for `M` steps starting from `s0`, using parameters in `params`.
Return value: named tuple with fields
    - `sequences`: alignment (or vector) of sequences
    - `tvals`: vector with the number of steps at each sample
    - `info`: information about the run
    - `params`: parameters of the run.


*Note*: this function is not very efficient if `M` is small.

Sampling details are determined by `parameters`, see `?SamplingParameters`.
Whether to use the genetic code is determined by the type of the `init` sequence:
it is used if `init::CodonSequence`, otherwise not.
"""
function mcmc_sample end

function mcmc_sample(
    g::PottsGraph,
    M::Int,
    s0::AbstractSequence,
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    verbose=0,
    progress_meter=true,
    alignment_output=true,
    translate_output=false,
)
    if !alignment_output && translate_output
        error("I have to implement this case")
    end

    @unpack Teq, burnin = params
    @argcheck Teq > 0 "Number of steps between samples `Teq` should be >0. Instead $(Teq)"
    @argcheck M > 0 "Number of samples `M` must be >0. Instead $M"
    tmp_check_alphabet_consistency(g, s0)
    verbose > 0 && @info """
          Sampling $M sequences using the following settings:
          - Type of sequence = $(supertype(s0))
          - Steps between samples = $(Teq)
          - burnin = $(burnin)
          - Step style = $(params.step_type)
          - Branch length style = $(params.step_meaning)
          - fraction of gap steps (if codon) = $(params.fraction_gap_step)
      """
    verbose > 0 && @info "Initial sequence: $s0"

    # vector for return value
    conf = copy(s0)
    S = similar([conf], M) # the sample
    tvals = Vector{Int}(undef, M) # number of steps at each sample

    # Holder if gibbs steps are used
    # the size of the holder depends on the symbols of the sequence: amino acids or codons
    gibbs_holder = get_gibbs_holder(s0)

    # Burnin
    verbose > 0 && println("Initializing with $(burnin) burnin iterations... ")
    burnin > 0 && mcmc_steps!(conf, g, burnin, params; rng, gibbs_holder)
    S[1] = copy(conf)
    tvals[1] = burnin

    # Sampling
    log_info = []
    progress = Progress(
        M - 1;
        barglyphs=BarGlyphs("[=> ]"),
        desc="Sampling: ",
        showspeed=true,
        enabled=progress_meter,
        dt=2,
    )
    time = @elapsed for m in 2:M
        # doing Teq steps on the current configuration
        _, proposed, performed = mcmc_steps!(
            conf, g, Teq, params; rng, gibbs_holder, verbose
        )
        # storing the result in S
        S[m] = copy(conf)
        tvals[m] = tvals[m-1] + Teq
        # misc.
        push!(
            log_info, (proposed=proposed, performed=performed, ratio=performed/proposed)
        )
        next!(progress; showvalues=[("steps", m + 1), ("total", M)])
    end
    verbose > 0 && @info "Sampling done in $time seconds"

    sequences = fmt_output(S, alignment_output, translate_output; names=tvals)
    return (; sequences, tvals, info=log_info, params=return_params(params, s0))
end
"""
    mcmc_sample(
        g::PottsGraph, M::Integer, params::SamplingParameters; init=:random_num, kwargs...)
    )

Sample `g` for `M` steps, using parameters in `params`.
Choose the initial sequence based on `init` and `g`.
See `PottsEvolver.get_init_sequence` for more information.
"""
function mcmc_sample(
    g::PottsGraph, M::Integer, params::SamplingParameters;
    init=:random_num, verbose=0, kwargs...
)
    s0 = get_init_sequence(init, g; verbose)
    return mcmc_sample(g, M, s0, params; verbose, kwargs...)
end

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
    sequence::AbstractSequence, g::PottsGraph, num_steps::Integer, p::SamplingParameters;
    rng=Random.GLOBAL_RNG, gibbs_holder=get_gibbs_holder(sequence), verbose=false,
)
    @argcheck length(sequence) > 0

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
        proposed >= max_tries && @warn """
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
    return if isgap(codon_ref, codon_alphabet)
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

        _metropolis_step!(s, ΔE, i, codon_new, aa_new, aa_ref)
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

            _metropolis_step!(s, ΔE, i, codon_new, aa_new, aa_ref)
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

#======================================================#
################### Initial sequence ###################
#======================================================#

"""
    get_init_sequence(s0, g::PottsGraph; kwargs...)

Try to guess a reasonable init sequence from `s0`:
- if `s0::AbstractSequence`, use a **copy** of it;
- if `s0::Symbol`, then it should be among `[:random_codon, :random_aa, :random_num]`;
  a random sequence of the corresponding type is created, using the length of `g`;
- if `s0` is a vector of integers, convert it to `AASequence`, `CodonSequence` or `NumSequence`;
  the conversion type depends on the alphabet of `g` and on the maximum element of `s0`.
"""
function get_init_sequence(s0::Symbol, g::PottsGraph; kwargs...)
    (; L, q) = size(g)
    return if s0 == :random_codon
        @argcheck q == length(aa_alphabet) """
            For sampling from `CodonSequence`, graph alphabet size must be $(length(aa_alphabet)).
            Instead $q.
            """
        CodonSequence(L)
    elseif s0 == :random_aa
        @argcheck q == length(aa_alphabet) """
            For sampling from `AASequence`, graph alphabet size must be $(length(aa_alphabet)).
            Instead $q.
            """
        AASequence(L)
    elseif s0 == :random_num
        NumSequence(L, q)
    else
        error(
            "Invalid symbol `init = $s0`. Options: `[:random_codon, :random_aa, :random_num]`",
        )
    end
end
get_init_sequence(s0::AbstractSequence, g; kwargs...) = copy(s0)
function get_init_sequence(s0::AbstractVector{<:Integer}, g; verbose=true)
    return if g.alphabet == aa_alphabet
        if maximum(s0) <= 21
            AASequence(s0)
        elseif 21 < maximum(s0) <= 65
            CodonSequence(s0)
        else
            error("Sequence $s0 incompatible with graph of size $(size(g))")
        end
    else
        NumSequence(s0)
    end
end

#=====================#
######## Utils ########
#=====================#

function fmt_output(
    sequences::AbstractVector{T}, alignment, translate;
    names=nothing, dict=false,
) where T<:CodonSequence
    return if alignment
        A = Alignment(sequences; names)
        translate ? genetic_code(A) : A
    elseif dict
        Dict{String, T}(name => seq for (name, seq) in zip(names, sequences))
    else
        sequences
    end
end
function fmt_output(
    sequences::AbstractVector{T}, alignment, translate;
    names=nothing, dict=false,
) where T<:AbstractSequence
    return if alignment
        Alignment(sequences; names)
    elseif dict
        Dict{String, T}(name => seq for (name, seq) in zip(names, sequences))
    else
        sequences
    end
end

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

get_gibbs_holder(::CodonSequence, T=FloatType) = zeros(T, 4)
get_gibbs_holder(::AASequence, T=FloatType) = zeros(T, 21)
get_gibbs_holder(s::NumSequence, T=FloatType) = zeros(T, s.q)

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

function tmp_check_alphabet_consistency(g::PottsGraph, s0::CodonSequence)
    if symbols(g.alphabet) != symbols(aa_alphabet)
        @warn """
            For now, sampling is only possible for graphs with the default alphabet $(aa_alphabet)
            Instead $(g.alphabet)
            """
    end
    return false
end
function tmp_check_alphabet_consistency(g::PottsGraph, s0::AASequence)
    if symbols(g.alphabet) != symbols(aa_alphabet)
        @warn """
            For now, sampling is only possible for graphs with the default alphabet $(aa_alphabet)
            Instead $(g.alphabet)
            """
    end
    return false
end
tmp_check_alphabet_consistency(g::PottsGraph, s0::AbstractSequence) = true

function return_params(p::SamplingParameters, s::T) where {T<:AbstractSequence}
    d = Dict()
    for field in propertynames(p)
        d[field] = getproperty(p, field)
    end
    d[:sequence_type] = T
    return d
end
