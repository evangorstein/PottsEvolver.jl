# just a wrapper to have <:TreeNodeData
@kwdef struct Sequence{T<:AbstractSequence} <: TreeNodeData
    seq::T = T(0)
end

Base.copy(S::Sequence) = Sequence(copy(S.seq))

function mcmc_sample_tree!(
    g::PottsGraph, tree::Tree, rootseq::AbstractSequence, params; kwargs...
)
    prepare_tree!(tree, rootseq)
    return mcmc_sample_tree!(g, tree, params; kwargs...)
end

"""
    mcmc_sample_single!(
        g::PottsGraph, tree::Tree{<:Sequence{S}}, params::SamplingParameters;
        rng=Random.GLOBAL_RNG, verbose=0,
    ) where S <: AbstractSequence

Sample one sequence per node of `tree`.
Expects that `tree` is already "ready": it has a non empty root sequence, and branch
lengths are integers.
"""
function mcmc_sample_tree!(
    g::PottsGraph,
    tree::Tree{<:Sequence{S}},
    params::SamplingParameters;
    rng=Random.GLOBAL_RNG,
    verbose=0,
) where {S<:AbstractSequence}
    # some settings
    rootseq = data(root(tree)).seq
    gibbs_holder = get_gibbs_holder(rootseq)

    # checks
    tmp_check_alphabet_consistency(g, rootseq)

    @argcheck all(n -> is_approx_integer(branch_length(n)), nodes(tree; skiproot=true)) """
        Branches of tree should be integers.\
        Instead $(map(branch_length, nodes(tree; skiproot=true)))
    """

    sample_children!(tree.root, g, params, gibbs_holder; rng, verbose)

    return tree
end

function sample_children!(node::TreeNode{<:Sequence}, g, params, gibbs_holder; kwargs...)
    for c in children(node)
        # copy the ancestral sequence
        s0 = copy(data(node).seq)
        # mcmc using it as an init
        nsteps = branch_length(c)
        mcmc_steps!(s0, g, nsteps, params; gibbs_holder, kwargs...) # from sampling.jl
        # copy result to child
        data!(node, Sequence(s0))
        # recursive
        sample_children!(c, g, params, gibbs_holder; kwargs...)
    end
    return nothing
end

#=====================#
######## Utils ########
#=====================#


is_approx_integer(x) = isapprox(x, round(x), rtol=1e-8)

function prepare_tree(tree::Tree, rootseq::S) where S <: AbstractSequence
    # convert to right type -- this makes a copy
    tree_copy = convert(Tree{Sequence{S}}, tree)
    # set root sequence
    data!(root(tree_copy), Sequence(copy(rootseq)))
    #
    round_tree_branch_length!(tree_copy)
    return tree_copy
end

function round_tree_branch_length!(tree::Tree; rtol=1e-3)
    for node in nodes(tree; skiproot=true)
        x = branch_length(node)
        if ismissing(x) || !isapprox(x, round(Int, x); rtol)
            error("Expected tree with (approximately) integer branch lengths. Instead $x")
        else
            branch_length!(node, round(Int, x))
        end
    end
    return tree
end

