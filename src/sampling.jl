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
    function SamplingParameters(step_type, branch_length_type, Teq, burnin)
        @assert branch_length_type in VALID_BRANCH_LENGTH_TYPES """
                `branch_length_type` should be in $VALID_BRANCH_LENGTH_TYPES.
                Instead $(branch_length_type).
            """
        @assert step_type in VALID_STEP_TYPES """
                `step_type` should be in $VALID_STEP_TYPES.
                Instead $(step_type).
            """
        return new(step_type, branch_length_type, Teq, burnin)
    end
end


"""
    mcmc_sample(
        g::PottsGraph, M::Integer, p::SamplingParameters;
        init = CodonSequence(size(g).L; source=:aa), # random amino acid sequence
        verbose=false
        rng = Random.GLOBAL_RNG,
    )
    mcmc_sample(g, M; init, verbose, rng, Teq = L, kwargs...)

Sample `g` for `M` steps starting from `init`.
Return a vector of sequences with the type of `init`, and a vector with the corresponding
    number of steps.

Sampling details are determined by `parameters`, see `?SamplingParameters`.
Whether to use the genetic code is determined by the type of the `init` sequence:
- it is used if `init::CodonSequence`
- it is not if `init::AASequence`

Second form: `kwargs` are passed to `SamplingParameters`.
"""
function mcmc_sample end

function mcmc_sample(
    g::PottsGraph, M::Integer, p::SamplingParameters;
    init = CodonSequence(size(g).L; source=:aa), # random amino acid sequence
    rng = Random.GLOBAL_RNG,
    verbose=false,
)
    verbose && println("""
        Sampling $M sequences using the following settings:
        - `Teq` = $(p.Teq)
        - Step style = $(p.step_type)
        - Branch length style = $(p.branch_length_type)
    """)

    L, q = size(g)
    conf = copy(init)
    S = similar([conf], M) # the sample
    tvals = Vector{Int}(undef, M) # the time values

    # Holder if gibbs steps are used
    # the size of the holder depends on the symbols of the sequence: amino acids or codons
    gibbs_holder = get_gibbs_holder(init)

    # Burnin
    verbose && println("Initializing with $(p.burnin) burnin iterations... ")
    for _ in 1:p.burnin
        mcmc_step!(conf, g, p; rng, gibbs_holder)
    end
    S[1] = copy(conf)
    tvals[1] = p.burnin
    # Sampling
    for m in 2:M
        for t in 1:p.Teq
            mcmc_step!(conf, g, p; rng, gibbs_holder)
        end
        S[m] = copy(conf)
        tvals[m] = tvals[m-1] + p.Teq
    end

    return S, tvals
end

function mcmc_sample(
    g::PottsGraph, M::Integer;
    init = CodonSequence(size(g).L; source=:aa), # random amino acid sequence
    rng = Random.GLOBAL_RNG,
    verbose=false,
    Teq = size(g).L,
    kwargs...
)
    return mcmc_sample(g, M, SamplingParameters(; Teq, kwargs...); init, rng, verbose)
end
#===================================================#
################## Codon sequences ##################
#===================================================#

"""
    mcmc_step!(s, g, p::SamplingParameters; kwargs...)

Determine the step to perform using `p`, then execute it.
"""
function mcmc_step!(
    s::AbstractSequence, g::PottsGraph, p::SamplingParameters;
    rng = Random.GLOBAL_RNG, gibbs_holder = get_gibbs_holder(s),
)
    step_func! = if p.step_type == :gibbs
        gibbs_step!
    elseif p.step_type == :metropolis
        metropolis_step!
    else
        error("Unknown `step_type` $(p.step_type)")
    end

    if p.branch_length_type == :proposed
        step_func!(s, g, p; rng, gibbs_holder)
    else
        max_tries = p.Teq * 10
        k = 0
        accepted = false
        while !accepted && k < max_tries
            accepted = step_func!(s, g, p; rng, gibbs_holder)[end]
            k += 1
        end
        k >= max_tries && @warn "$max_tries steps attempted without acceptance. Giving up."
    end

    return s
end

#===============================================================#
###################### CodonSequence steps ######################
#===============================================================#

function gibbs_step!(s::CodonSequence, g::PottsGraph, p::SamplingParameters; kwargs...)
    # need to implement gap step here, with some parameter coming from p probably
    return aa_gibbs_step!(s, g; kwargs...)
end

function metropolis_step!(s::CodonSequence, g::PottsGraph, p::SamplingParameters; kwargs...)
    error("Not implemented yet")
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
            # ΔE is _minus_ the energy --> softmax(ΔE) is the Boltzmann distribution
            ΔE = g.h[aanew,i] - g.h[aaref,i]
            for j in 1:length(s)
                if j != i
                    ΔE += g.J[aanew, s.aaseq[j], i, j] - g.J[aaref, s.aaseq[j], i, j]
                end
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
    change = (new_aas[c] != aaref)
    return i, b, s[i], change
end

function gap_gibbs_step!(s::CodonSequence, g::PottsGraph; rng = Random.GLOBAL_RNG)
    error("Not implemented")
    return false
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
