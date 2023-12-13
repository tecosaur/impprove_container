mutable struct ConsensusModel{M, S} <: MMI.Model where {M <: MMI.Model, S <: ResamplingStrategy}
    base::M
    resampler::S
    n::Int
    oob::Bool
    importance::Bool
    thread::Bool
    distribute::Bool
end

ConsensusModel(base::M, resampler::ResamplingStrategy,
               n::Int=100; oob::Bool=false, importance::Bool=false,
               thread::Bool=false, distribute::Bool=false) where { M <: MMI.Model } =
    ConsensusModel(base, resampler, n, oob, importance, thread, distribute)

using Base.GC

function MMI.fit(cm::ConsensusModel, verbosity::Int, X, y)
    p = Progress(cm.n, enabled=verbosity > 0)
    X1, y1 = MMI.reformat(cm.base, X, y)
    submachines =
        (train_test_pairs(cm.resampler, axes(X1, 1)) for _ in 1:cm.n) |>
        Map(pairs -> if length(pairs) == 1
                first(pairs)
            else
                Iterators.map(Iterators.flatten, zip(pairs...))
            end) |>
        Map(function ((train, test),)
                trainX, trainY = @views X1[train, :], y1[train]
                mach = machine(cm.base, trainX, trainY)
                fit!(mach, verbosity=verbosity-1)
                oob, testX, testY = if cm.oob
                    stest = sort(test)
                    testX = @view X1[stest, :]
                    (stest, predict(mach, testX)), testX, y1[test]
                else
                    nothing, trainX, trainY
                end
                importances = if cm.importance
                    consensus_feature_importances(mach, testX, testY)::Vector{<:Pair}
                end
                next!(p)
                (; fitres=mach.fitresult, test, importances)
            end) |>
                if cm.distribute dcollect
                elseif cm.thread tcollect
                else collect end
    oob = if cm.oob
        merge(aggregated_oob(cm, verbosity, getproperty.(submachines, :fitres),
                             X, getproperty.(submachines, :test)),
              (; actual=MLJBase.CategoricalArrays.unwrap.(y)))
    end
    fitresult = (; components=getproperty.(submachines, :fitres),
                 # testsets=getproperty.(submachines, :test),
                 # weights = ones(Float64, length(components)),
                 importances = Dict(
                     sort(first.(first(submachines).importances)[:]) .=>
                         (mapreduce(subm -> last.(sort(subm.importances[:], by=first)),
                                    +, submachines) ./ length(submachines))),
                 oob)
    (fitresult, nothing, report)
end

# If `X` is a dataframe, this tends to cause type instability issues,
# so it's worth converting it to a MatrixTable.
MMI.fit(cm::ConsensusModel, verbosity::Int, X::DataFrame, y) =
    MMI.fit(cm, verbosity, MatrixTable(X), y)

"""
    consensus_feature_importances(mach::Machine, testX, testY)
Obtain the feature importances for `mach`, optionally baesd on `testX` and `testY`.
"""
consensus_feature_importances(mach::Machine, _, _) =
    MMI.feature_importances(mach)

function aggregated_oob(cm::ConsensusModel, verbosity::Int, cfitres::Vector, X, testsets::Vector)
    verbosity >= 1 && @info "Aggregating OOB predictions"
    oob = [MMI.predict(cm.base, fitres, X[test, :]) for (fitres, test) in zip(cfitres, testsets)]
    results = missings(typeof(first(first(oob))), size(X, 1), length(oob))
    for (col, (rows, values)) in enumerate(zip(testsets, oob))
        results[rows, col] = values
    end
    eachrow(results) |> Map(skipmissing) |> Map(collect) |>
        Map(vals -> if length(vals) > 0 mean(vals) else missing end) |> tcollect
end

function MMI.predict(cm::ConsensusModel, fitresult, new)
    mean.(fitresult.components |> Map(c -> MMI.predict(cm.base, c, new)) |> tcollect)
end
