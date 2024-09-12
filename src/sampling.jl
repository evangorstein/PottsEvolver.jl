const VALID_STEP_TYPES = [:gibbs, :metropolis]
const VALID_BRANCH_LENGTH_TYPES = [:accepted, :proposed]

"""
    mutable struct SamplingParameters

Construct using keyword arguments:
```
step_type::Symbol = :gibbs
branch_length::Symbol = :accepted
Teq::Int
burnin::Int = 5*Teq
```

- `Teq` is measured in swaps: attempted (or accepted) change of one sequence position.
- `burnin`: number of steps starting from the initial sequence before the first `Teq` are
    made.
- if `branch_length_type` is `:accepted`, only accepted MCMC steps will count towards
    equilibration. If it is `:proposed`, all steps count.
    *e.g.*: `Teq = 10`, starting from sequence `s`. If `:accepted` is used, then the next
    sample is taken after `10` accepted MCMC steps take place. If `:proposed`, `10` calls
    to the stepping function are enough.
    **Important**: only steps that lead to amino acid replacements are considered as
    accepted.
"""
@kwdef mutable struct SamplingParameters
    step_type::Symbol = :gibbs
    branch_length_type::Symbol = :accepted
    Teq::Int
    burnin::Int = 5*Teq
    fraction_gap_step = 0.9
    function SamplingParameters(step_type, branch_length_type, Teq, burnin, fraction_gap_step)
        @assert branch_length_type in VALID_BRANCH_LENGTH_TYPES """
                `branch_length_type` should be in $VALID_BRANCH_LENGTH_TYPES.
                Instead $(branch_length_type).
            """
        @assert step_type in VALID_STEP_TYPES """
                `step_type` should be in $VALID_STEP_TYPES.
                Instead $(step_type).
            """
        return new(step_type, branch_length_type, Teq, burnin, fraction_gap_step)
    end
end


"""
    mcmc_sample(
        g::PottsGraph, M::Integer, p::SamplingParameters;
        init = :random,
        verbose=false
        rng = Random.GLOBAL_RNG,
    )
    mcmc_sample(g, M; init, verbose, rng, Teq = L, kwargs...)

Sample `g` for `M` steps starting from `init`.
Return a vector of sequences with the type of `init`, and a vector with the corresponding
    number of steps.

*Note*: this function is not very efficient if `M` is small.

Sampling details are determined by `parameters`, see `?SamplingParameters`.
Whether to use the genetic code is determined by the type of the `init` sequence:
- it is used if `init::CodonSequence`
- it is not if `init::AASequence` or `init::NumSequence`

If the `init=:random` argument is used, the initial sequence will be
- a random `CodonSequence` if `size(g).q == 21`, by a call to `CodonSequence(sequence_length)`
- a `NumSequence` otherwise

Second form: `Teq` and `kwargs` are passed to `SamplingParameters`.
"""
function mcmc_sample end

function mcmc_sample(
    g::PottsGraph, M::Integer, p::SamplingParameters;
    init = :random,
    rng = Random.GLOBAL_RNG,
    verbose=false,
    progress_meter = true,
)
    @assert p.Teq > 0 "Number of steps between samples `Teq` should be >0. Instead $(p.Teq)"
    verbose && @info """
        Sampling $M sequences using the following settings:
        - `Teq` = $(p.Teq)
        - Step style = $(p.step_type)
        - Branch length style = $(p.branch_length_type)
    """

    L, q = size(g)

    # initial sequence and vector for return value
    s0 = (init == :random) ? random_init_sequence(g, nothing) : init
    conf = copy(s0)
    S = similar([conf], M) # the sample
    tvals = Vector{Int}(undef, M) # the time values

    # Holder if gibbs steps are used
    # the size of the holder depends on the symbols of the sequence: amino acids or codons
    gibbs_holder = get_gibbs_holder(s0)

    # Burnin
    verbose && println("Initializing with $(p.burnin) burnin iterations... ")
    p.burnin > 0 && mcmc_steps!(conf, g, p.burnin, p; rng, gibbs_holder)
    S[1] = copy(conf)
    tvals[1] = p.burnin

    # Sampling
    log_info = []
    progress = Progress(
        M-1;
        barglyphs = BarGlyphs("[=> ]"),
        desc="Sampling: ",
        showspeed=true,
        enabled = progress_meter,
    )
    time = @elapsed for m in 2:M
        # doing p.Teq steps on the current configuration
        _, proposed, performed = mcmc_steps!(conf, g, p.Teq, p; rng, gibbs_holder, verbose)
        # storing the result in S
        S[m] = copy(conf)
        tvals[m] = tvals[m-1] + p.Teq
        # misc.
        push!(log_info, (proposed=proposed, performed=performed, ratio=performed/proposed))
        next!(
            progress;
            showvalues = [("steps", m+1), ("total", M)],
        )
    end
    verbose && @info "Sampling done in $time seconds"

    return Alignment(S; alphabet=g.alphabet, names = tvals), log_info
end

function mcmc_sample(
    g::PottsGraph, M::Integer;
    init = :random, # random amino acid sequence
    rng = Random.GLOBAL_RNG,
    verbose = false,
    progress_meter = true,
    Teq = size(g).L,
    kwargs...
)
    return mcmc_sample(g, M, SamplingParameters(; Teq, kwargs...); init, rng, verbose)
end

"""
    mcmc_steps!(
        s::AbstractSequence, g, num_steps, p::SamplingParameters;
        gibbs_holder, kwargs...
    )

Perform `num_steps` MCMC steps starting from sequence `s` and using graph `g`.
The step type (`:gibbs`, `:metropolis`) and the interpretation of `num_steps`
    (`:accepted`, `:proposed`) is set using `p` (see `?SamplingParameters`).
Modifies the input sequence `s` and returns it.
"""
function mcmc_steps!(
    s::AbstractSequence, g::PottsGraph, num_steps::Integer, p::SamplingParameters;
    rng = Random.GLOBAL_RNG, gibbs_holder = get_gibbs_holder(s), verbose=false,
)
    @debug "Performing $num_steps ($(p.branch_length_type)) steps of type $(p.step_type)"

    step_func! = if p.step_type == :gibbs
        gibbs_step!
    elseif p.step_type == :metropolis
        metropolis_step!
    else
        error("Unknown `step_type` $(p.step_type)")
    end

    proposed = 0
    performed = 0
    if p.branch_length_type == :proposed
        for _ in 1:num_steps
            step_func!(s, g, p; rng, gibbs_holder)
            proposed += 1
            performed += 1
        end
    else
        min_acceptance_rate = 0.01
        max_tries = num_steps / min_acceptance_rate
        while performed < num_steps && proposed < max_tries
            step_result = step_func!(s, g, p; rng, gibbs_holder)
            if p.branch_length_type == :accepted && step_result.accepted
                performed += 1
            elseif p.branch_length_type == :changed && step_results.changed
                performed += 1
            end
            proposed += 1
        end
        proposed >= max_tries && @warn """
            $max_tries steps attempted with only $performed performed. Giving up.
            """
    end

    return s, proposed, performed
end

#=
MCMC step functions
the return value should be a named tuple of the form
`(i, new_state, accepted, changed)`, plus extra info if needed
`accepted` says whether the step was accepted (always true for gibbs)
`change` is specific to codon sequences, and says whether the amino acid state was changed.
=#


#===============================================================#
###################### CodonSequence steps ######################
#===============================================================#

function gibbs_step!(s::CodonSequence, g::PottsGraph, params::SamplingParameters; kwargs...)
    p = params.fraction_gap_step
    return if rand() < p
        gap_metropolis_step!(s, g; kwargs...)
    else
        aa_gibbs_step!(s, g; kwargs...)
    end
end

function metropolis_step!(s::CodonSequence, g::PottsGraph, p::SamplingParameters; kwargs...)
    error("Not implemented")
end

"""
    aa_gibbs_step!(s::CodonSequence, g::PottsGraph; kwargs...)

Change one coding codon in `s` into another coding codon.
"""
function aa_gibbs_step!(
    s::CodonSequence, g::PottsGraph;
    rng = Random.GLOBAL_RNG, gibbs_holder = zeros(Float64, 4),
)
    p = gibbs_holder

    # new_codons and new_aas are arrays but should _NOT_ be mutated
    i, b, new_codons, new_aas = pick_aa_mutation(s; rng)
    @debug "New possible codons / aas" new_codons new_aas
    aaref = s.aaseq[i]
    for (a, aanew) in enumerate(new_aas)
        if aanew == aaref
            p[a] = 0. # energy for now
        else
            # ΔE is minus the energy --> softmax(ΔE) is the Boltzmann distribution
            # ~ high ΔE is more probable
            ΔE = - aa_degeneracy(aanew) + aa_degeneracy(aaref) # log-degeneracy of gen code
            ΔE += g.h[aanew,i] - g.h[aaref,i]
            for j in Iterators.filter(!=(i), 1:length(s))
                ΔE += g.J[aanew, s.aaseq[j], i, j] - g.J[aaref, s.aaseq[j], i, j]
            end
            p[a] = g.β*ΔE
        end
    end
    if length(new_aas) < length(p)
        p[(length(new_aas)+1):end] .= -Inf
    end
    softmax!(p)
    @debug "Gibbs conditional probability" p
    c = wsample(p)
    s[i] = new_codons[c]
    s.aaseq[i] = new_aas[c]
    changed = (new_aas[c] != aaref)
    accepted = true
    return (i=i, b=b, new_state=s[i], accepted=accepted, changed=changed)
end

function gap_metropolis_step!(
    s::CodonSequence, g::PottsGraph;
    rng = Random.GLOBAL_RNG, kwargs...,
)
    i = rand(1:length(s))
    codon_ref = s.seq[i]
    aa_ref = s.aaseq[i]
    # di Bari et. al.
    β = 1/n_aa_codons # aa to gap transition
    return if isgap(codon_ref)
        # Pick any non gap codon and try to mutate to it
        while (codon_new = rand(codon_alphabet.char_to_index)[2]; !iscoding(codon_new)) end
        aa_new = genetic_code(codon_new)
        ΔE = - aa_degeneracy(aa_new) + aa_degeneracy(aa_ref)
        ΔE += g.h[aa_new, i] - g.h[aa_ref, i]
        for j in Iterators.filter(!=(i), 1:length(s))
            ΔE += g.J[aa_new, s.aaseq[j], i, j] - g.J[aa_ref, s.aaseq[j], i, j]
        end

        _metropolis_step!(s, ΔE, i, codon_new, aa_new, aa_ref)
    else
        # Try to replace s[i] by the gap codon
        if rand() < 1-β
            # with probability 1-β do nothing
            # return value below should match _metropolis_step!(::CodonSequence, ...)
            (i=i, new_state=s[i], accepted=false, changed=false)
        else
            # Metropolis move: replace s[i] with gap codon
            codon_new = codon_alphabet(gap_codon) # gap_codon defined in codons.jl
            aa_new = genetic_code(codon_new)
            ΔE = aa_degeneracy(aa_ref) # degeneracy of gap state (aa_new) is 0
            ΔE += g.h[aa_new, i] - g.h[aa_ref, i]
            for j in Iterators.filter(!=(i), 1:length(s))
                ΔE += g.J[aa_new, s.aaseq[j], i, j] - g.J[aa_ref, s.aaseq[j], i, j]
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
        (i=i, new_state=s[i], accepted=true, changed=aa_new==aa_ref)
    else
        (i=i, new_state=s[i], accepted=false, changed=false)
    end
end

#===================================================#
################## Non-codon steps ##################
#===================================================#

function gibbs_step!(
    s::AbstractSequence, g::PottsGraph, p::SamplingParameters;
    rng = Random.GLOBAL_RNG, gibbs_holder = get_gibbs_holder(S)
)
    p = gibbs_holder
    q = length(gibbs_holder)

    i = rand(1:length(s))
    a_ref = s[i]
    for a in 1:q
        if a == a_ref
            p[a] = 0.
        else
            ΔE = g.h[a,i] - g.h[a_ref,i]
            for j in 1:length(s)
                if j != i
                    ΔE += g.J[a, s[j], i, j] - g.J[a_ref, s[j], i, j]
                end
            end
            p[a] = g.β*ΔE
        end
    end

    softmax!(p)
    @debug "Gibbs conditional probability" p
    c = wsample(p)
    s[i] = wsample(p)
    change = (a_ref != s[i])
    return i, s[i], change
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
function pick_aa_mutation(s::CodonSequence; rng = Random.GLOBAL_RNG)
    max_cnt = 1000
    cnt = 1
    # Pick i and b, and see if we can mutate the codon (if it's a gap we can't)
    @label again
    cnt += 1
    i = rand(1:length(s))
    b = rand(1:3)
    codon = s[i]
    new_codons, new_aas = accessible_codons(codon, b)
    if cnt < max_cnt && (isnothing(new_codons) || isnothing(new_aas))
        @goto again
    else
        return i, b, new_codons, new_aas
    end
end

get_gibbs_holder(::CodonSequence) = zeros(Float64, 4)
get_gibbs_holder(::AASequence) = zeros(Float64, 21)
get_gibbs_holder(s::NumSequence) = zeros(Float64, s.q)

function softmax!(X)
    # thanks chatGPT!
    # Compute the maximum value in X to ensure numerical stability
    max_val = maximum(X)
    Z = 0.
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

"""
    random_init_sequence(g::PottsGraph)

Construct a random sequence of the type:
- `CodonSequence` if `q == 21` (call to `CodonSequence(sequence_length)`)
- `NumSequence` otherwise
"""
function random_init_sequence(g::PottsGraph, p)
    # p not used, might use it in the future so I'm leaving it here
    L, q = size(g)
    return q == 21 ? CodonSequence(L) : NumSequence(L, q)
end
