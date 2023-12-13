using LinearAlgebra
using EvalMetrics

function confusion_matricies(targets, predictions)
    thresholds = predictions |> copy
    sort!(thresholds)
    unique!(thresholds)
    (EvalMetrics.ConfusionMatrix(targets, predictions, thresholds),
     thresholds)
end

function πinterpolated(πs, αs, α)
    p = sortperm(αs)
    left = something(findlast(<=(α), αs[p]), 1)
    right = findfirst(>(α), αs[p])
    πs[p[left]] + (α - αs[p[left]]) *
        (πs[p[right]] - πs[p[left]]) / (αs[p[right]] - αs[p[left]]),
    p[left]
end

function scorerates(targets::AbstractVector{Bool}, scores::Vector{Float64})
    confmats, thresholds = confusion_matricies(targets, scores)
    rates = Dict(
        :thresholds => thresholds,
        :p => getproperty.(confmats, :p),
        :n => getproperty.(confmats, :n),
        :tp => getproperty.(confmats, :tp),
        :tn => getproperty.(confmats, :tn),
        :fp => getproperty.(confmats, :fp),
        :fn => getproperty.(confmats, :fn))
    rates[:precision] = @. rates[:tp] / (rates[:tp] + rates[:fp])
    rates[:recall] = @. rates[:tp] / rates[:p]
    rates[:α] = @. rates[:fp] / rates[:n]
    rates[:π] = @. rates[:tp] / rates[:p]
    rates[:f1] = @. 2 * rates[:precision] * rates[:recall] /
                    (rates[:precision] + rates[:recall])
    rates
end

"""
    scoreproba(targets::AbstractVector{Bool}, scores::Vector{Float64})
Calculate various metrics which are indicative of how well `scores` predict
`targets`, namely:
- `auROC`, the area under the reciver-operator curve
- `auPRC`, the area under the precision-recall curve
- `auPRCΔtpr`, the auPRC less the true positive rate
- `f1max`, the maximum attainable f1 value
- `topNtp`, the number of the top `N` scored entries that are
  true positives (for `N` = 20, 50, 100, 200)
- `πα⁻¹N`, the true positive rate (π) when the false positive rate (α)
  is at most 1 in `N` (for `N` = 100, 200, 500, 1000, 2000, 5000, 10000, 50000)
- `truepositive`, the total number of true positives
- `positiveratio`, the proportion of true positives (TPR)
- `nvals`, the number of values
"""
function scoreproba(targets::AbstractVector{Bool}, scores::Vector{Float64})
    results = Dict{Symbol, Any}()
    rates = scorerates(targets, scores)

    results[:auROC] = auc_trapezoidal(rates[:α], rates[:π])
    results[:auPRC] = auc_trapezoidal(rates[:recall], rates[:precision])
    results[:auPRCΔtpr] = results[:auPRC] - mean(targets)
    results[:f1max] = maximum(filter(!isnan, rates[:f1]))
    results[:truepositive] = sum(targets)
    results[:nvals] = length(targets)
    results[:positiveratio] = mean(targets)

    targetbyscore = targets[sortperm(scores, rev=true)]
    for n in (20, 50, 100, 200)
        m = min(n, results[:truepositive])
        count = sum(targetbyscore[1:m])
        results[Symbol("top$(n)tp")] = round(Int, count * n/m)
    end

    for α⁻¹ in (100, 200, 500, 1000, 2000, 5000, 10000, 50000)
        results[Symbol("πα⁻¹$(α⁻¹)")], _ =
            πinterpolated(rates[:π], rates[:α], 1/α⁻¹)
    end

    results
end

"""
    πα⁻¹rescale(πα⁻¹::Number, p::Integer, n::Integer, m::Integer)
Rescale a `πα⁻¹` value from a total population size of `p`+`n` to `p`+`m`,
where `p` is the number of positives, `n` is the initial number of negatives,
and `m` is the rescaled number of negatives.

This relies on the assumption that `n` and `m` represent identically distributed
samples.

## Example

"""
function πα⁻¹rescale(πα⁻¹::Number, p::Integer, n::Integer, m::Integer)
    tp = πα⁻¹rescale * p
    println("TODO")
end

"""
    scoreproba(d::Dict)
Run `scoreproba(d[:predictions])`
"""
scoreproba(d::Dict) = scoreproba(d[:predictions])

"""
    scoreproba(predictions::DataFrame)
Run `scoreproba(predictions.pathogenic, predictions.score)`, excluding
NaN values.
"""
function scoreproba(predictions::DataFrame)
    nonnan_preds = filter(:score => !isnan, predictions)
    scoreproba(nonnan_preds.pathogenic, nonnan_preds.score)
end

function scoreproba(results::NamedTuple{(:predicted, :actual), <:Tuple{<:MLJBase.CategoricalDistributions.UnivariateFiniteVector, <:AbstractVector{Bool}}})
    scoreproba(results.actual, pdf.(results.predicted, true))
end
