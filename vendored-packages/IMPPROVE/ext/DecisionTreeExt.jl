module DecisionTreeExt

# We special-case DecisionTreeClassifier because
# the default implementation is attrociously inefficient.

using IMPPROVE

using MLJBase
using MLJModelInterface
if isdefined(Base, :get_extension)
    using MLJDecisionTreeInterface
else
    using ..MLJDecisionTreeInterface
end
using Base.Threads

const MMI = MLJModelInterface
const DT = MLJDecisionTreeInterface.DT

function IMPPROVE.aggregated_oob(::ConsensusModel{DecisionTreeClassifier}, verbosity::Int, cfitres::Vector, X, testsets::Vector)
    verbosity >= 1 && @info "Aggregating OOB predictions"
    classes_seen, integers_seen = first(cfitres)[2:3]
    Xmat = MMI.matrix(X)
    pred_list = [zeros(Float64, size(Xmat, 1), length(classes_seen))
                 for _ in 1:nthreads()]
    @threads for i in 1:length(cfitres)
        fitres, test = cfitres[i], testsets[i]
        pred_list[threadid()][test, :] .+= DT.apply_tree_proba(
            first(fitres), view(Xmat, test, :), integers_seen)
    end
    global stuff = (; cfitres, testsets, Xmat, pred_list)
    preds = reduce(+, pred_list)
    counts = sum(preds; dims=2)
    (; # predacc, predcounts,
     predicted = MMI.UnivariateFinite(classes_seen, preds ./ counts))
end

function IMPPROVE.consensus_feature_importances(mach::Machine{DecisionTreeClassifier}, testX, testY)
    sort(vec(mach.fitresult[4] .=>
        mean(tree_permutation_importance(
            mach.fitresult[1], Matrix(testX), MMI.int(testY), mach.fitresult[3],
            5, 1000), dims=2)), by=last)
end

function MMI.predict(::ConsensusModel{DecisionTreeClassifier}, fitresult, X)
    Xmat = MMI.matrix(X)
    cfitres = fitresult.components
    classes_seen, integers_seen = first(cfitres)[2:3]
    @assert all(Ref(classes_seen) .== getindex.(cfitres, 2))
    @assert all(Ref(integers_seen) .== getindex.(cfitres, 3))
    sum_scores =
        mapreduce(tree -> DT.apply_tree_proba(tree, Xmat, integers_seen),
                  +, first.(cfitres))
    MMI.UnivariateFinite(classes_seen, sum_scores ./ length(cfitres))
end

"""
    tree_permutation_importance(tree::DT.Root,
                                X::Matrix, y::Vector, levels::Vector, repeats::Int=5)

Permutation importance using a stratified brier score.
"""
function tree_permutation_importance(tree::DT.Root,
                                     X::AbstractMatrix, y::Vector, levels::Vector,
                                     repeats::Int=5, maxclassrep::Number=Inf)
    yclassindices = map(l -> findall(l .== y), levels)
    if any(length.(yclassindices) .> maxclassrep)
        yclassindices = map(yinds -> yinds[1:min(length(yinds), maxclassrep)],
                            yclassindices)
        allyinds = sort(Iterators.flatten(yclassindices) |> collect)
        yindmap = zeros(Int, length(y))
        yindmap[allyinds] .= 1:length(allyinds)
        yclassindices = map(yinds -> yindmap[yinds], yclassindices)
        X = X[allyinds, :]
    end
    function stratified_brier(ŷ)
        map(enumerate(yclassindices)) do (n, yinds)
            if isempty(yinds) 0.0 else
                sum((1 .- ŷ[yinds, n]).^2) / length(yinds)
            end
        end |> mean
    end
    score(X) = stratified_brier(DT.apply_tree_proba(tree, X, levels))
    featids(n::DT.Node) =
        vcat(n.featid, featids(n.left), featids(n.right))
    featids(::DT.Leaf) = Int64[]
    featids(r::DT.Root) = featids(r.node)
    baseline = score(X)
    features = featids(tree) |> unique
    results = zeros(Float64, size(X, 2), repeats)
    for r in 1:repeats
        for f in features
            f_col = copy(X[:, f])
            shuffle!(@view X[:, f])
            results[f, r] = max(0.0, baseline - score(X))
            X[:, f] = f_col
        end
    end
    if any(isnan, results)
        @info "! NaN importance"
        global permimp = (; tree, X, y, levels, repeats, maxclassrep)
    end
    results
end

end
