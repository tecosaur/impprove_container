using Random
using StatsBase

import StatsBase.sample

"""
A wrapper around some data of type `D`, which can be used
to produce samples of the data.
"""
abstract type AbstractSampler{D} end

"""
    traintest(as::AbstractSampler) -> Tuple(train, test)
Obtain a pair (train, test) of sample sets of `as`.
"""
function traintest end

"""
    traintest(as::AbstractSampler, n::Integer) -> Iterator -> Tuple(train, test)
Obtain an iterator providing `n` (train, test) pairs.
"""
traintest(as::AbstractSampler, n::Integer) =
    Iterators.map(_ -> traintest(as), 1:n)

"""
    traintest(n::Integer) -> Function -> Iterator -> Tuple(train, test)
Obtain a function that when provided called on an AbstractSampler will
produce an iterator providing `n` (train, test) pairs.
"""
traintest(n::Integer) = Base.Fix2(traintest, n)

"""
    sample(as::AbstractSampler) -> sample
Obtain a sample from `as`, by default the "train" pair from `traintest`.
"""
sample(as::AbstractSampler) =
    first(traintest(as))

"""
    sample(as::AbstractSampler, n::Integer) -> Iterator -> sample
Obtain an iterator providing `n` samples.
"""
sample(as::AbstractSampler, n::Integer) =
    Iterators.map(_ -> sample(as), 1:n)

"""
    sample(n::Integer) -> Function -> Iterator -> sample
Obtain a function that when provided called on an AbstractSampler will
produce an iterator providing `n` samples.
"""
sample(n::Integer) = Base.Fix2(sample, n)

# General purpose functions

features(s::AbstractSampler) = axes(s.data, 1)
features(s::AbstractSampler{<:AbstractDataFrame}) = names(s.data)
Base.length(s::AbstractSampler) = length(s.data)
Base.size(s::AbstractSampler) = size(s.data)
Base.size(s::AbstractSampler{<:GroupedDataFrame}) = size(s.data.parent)
Base.size(s::AbstractSampler{<:Vector{<:GroupedDataFrame}}) = size(first(s.data).parent)
Base.size(s::AbstractSampler, dim::Integer) = size(s)[dim]

Base.show(io::IO, a::AbstractSampler) =
    print(io, sprint(show, typeof(a)), size(a))

# ------------
# Split/processed sampling
# ------------

abstract type AbstractSplitSampler{D} <: AbstractSampler{D} end

struct SplitSampler{D, R} <: AbstractSplitSampler{D} where {
    R <: Union{Symbol, <:Integer, <:AbstractString,
               Vector{<:Union{Symbol, <:Integer, <:AbstractString}}}}
    sampler::AbstractSampler{D}
    explanatory::Vector{<:Union{<:Integer, Symbol, <:AbstractString}}
    response::R
end

features(s::SplitSampler) = features(s.sampler)
Base.length(s::SplitSampler) = length(s.sampler)
Base.size(s::SplitSampler) = size(s.sampler)

function splitsample(s::AbstractSampler{D}, response::Union{<:R, Vector{<:R}};
                     exclude::Union{<:E, Vector{<:E}}=Symbol[],
                     include::Vector{<:I}=names(
                         if s.data isa Vector first(s.data) else s.data end)) where {
                         R <: Union{Symbol, <:AbstractString},
                         E <: Union{Symbol, <:AbstractString},
                         I <: Union{Symbol, <:AbstractString},
                         D <: Union{<:AbstractDataFrame, <:GroupedDataFrame, <:Vector{<:GroupedDataFrame}}}
    if s isa SplitSampler
        throw(ArgumentError("Cannot split a SplitSampler"))
    else
        evars = setdiff(string.(vcat(include)), string.(vcat(response)), string.(exclude))
        SplitSampler{D, typeof(Symbol.(response))}(
            s, Symbol.(evars), Symbol.(response))
    end
end

function splitsample(s::SplitSampler{<:Union{<:AbstractDataFrame, <:GroupedDataFrame, <:Vector{<:GroupedDataFrame}}, R},
                     sample::AbstractDataFrame) where {R}
    view(sample, :, s.explanatory),
    if R isa Vector
        view(sample, :, s.response)
    else
        sample[!, s.response]
    end
end

sample(s::SplitSampler) = splitsample(s, sample(s.sampler))
traintest(s::SplitSampler) = splitsample.(Ref(s), traintest(s.sampler))

# ------------
# Bootstrap
# ------------

struct Bootstrap{D} <: AbstractSampler{D}
    data::D
    parameters::NamedTuple
end

Bootstrap(data::Any; kwargs...) =
    Bootstrap(data, NamedTuple(kwargs))

sample(b::Bootstrap{<:OrdinalRange}) = rand(b.data, length(b.data))
function traintest(b::Bootstrap{<:OrdinalRange})
    train = sample(b)
    train, setdiff(b.data, train)
end

sample(b::Bootstrap{<:Union{<:AbstractDataFrame, <:AbstractMatrix}}) =
    b.data[rand(axes(b.data, 1), size(b.data, 1)), :]
function traintest(b::Bootstrap{<:Union{<:AbstractDataFrame, <:AbstractMatrix}})
    train_rows = rand(axes(b.data, 1), size(b.data, 1))
    view(b.data, train_rows, :),
    view(b.data, setdiff(axes(b.data, 1), train_rows), :)
end

"""
    traintest(b::Bootstrap{GroupedDataFrame})
## Relevant parameters
- `grab`, which rows from each group should be grabbed, one of:
  - `:group` (default), the entire group
  - `:singular`, a single row from each group
  - `:samesingular`, a particular (always the same) single row from each group
- `traingrab`, a train-set specific version of `grab`
- `testgrab`, a train-set specific version of `grab`
- `truncate`, a maximum train group size beyond which data should be chopped off
- `truncatetotest`, whether truncated items from the train set should be moved
  into the test set.
- `indices` (Bool=false), whether to return the row indices instead of the rows
"""
function traintest(b::Bootstrap{<:GroupedDataFrame})
    train_groups = rand(axes(b.data, 1),
                        min(length(b.data),
                            get(b.parameters, :truncate, typemax(Int))))
    test_groups = setdiff(axes(b.data, 1), train_groups)
    if get(b.parameters, :truncate, nothing) !== nothing &&
        get(b.parameters, :truncatetovoid, false)
        nstrip = max(0, length(b.data) - b.parameters.truncate)
        rand_strip = rand(axes(b.data, 1), nstrip)
        test_groups = filter(r -> r ∉ rand_strip, test_groups)
    end
    train_grab = get(b.parameters, :traingrab,
        get(b.parameters, :grab, :group))
    test_grab = get(b.parameters, :testgrab,
        get(b.parameters, :grab, :group))
    group_rows = zip(b.data.starts, b.data.ends) .|> r -> b.data.idx[r[1]:r[2]]

    function grab_rows(groups, grab)
        if grab === :singular
            rand.(group_rows[groups])
        elseif grab === :samesingular
            rand.(group_rows)[groups]
        elseif grab === :group
            group_rows[groups] .|> collect |> Iterators.flatten |> collect
        else
            throw(ArgumentError("$b contains invalid grab parameter: $(grab)"))
        end
    end

    train_rows = grab_rows(train_groups, train_grab)
    test_rows = grab_rows(test_groups, test_grab)

    if get(b.parameters, :indices, false) === true
        train_rows, test_rows
    else
        view(b.data.parent, train_rows, :),
        view(b.data.parent, test_rows, :)
    end
end

"""
    traintest(b::Bootstrap{Vector{GroupedDataFrame}})
## Relevant parameters
- `scheme`, the sampling scheme to use, should be one of:
  - `:stratified` (default)
  - `:truncated`
- `grab`, which rows from each group should be grabbed, one of:
  - `:group` (default), the entire group
  - `:singular`, a single row from each group
  - `:samesingular`, a particular (always the same) single row from each group
- `traingrab`, a train-set specific version of `grab`
- `testgrab`, a train-set specific version of `grab`
- `truncationratio` (Number=1), the applied maximum train group size,
  relative to the minimum train group size
- `truncatetotest` (Bool=false), relocate truncated entries to the test set
- `ensurenonempty` (Bool=true), re-run if either the test or train group is empty
- `indices` (Bool=false), whether to return the row indices instead of the rows
"""
function traintest(b::Bootstrap{<:Vector{<:GroupedDataFrame}})
    function some(predicate::Function, sequence)
        for s in sequence
            r = predicate(s)
            if !isnothing(r)
                return r
            end
        end
    end
    if !all(gdf -> gdf.parent ===(first(b.data).parent), b.data)
        throw(ArgumentError("$b GroupedDataFrames do not all share the same parent"))
    end
    groups = filter(g -> length(g) > 0, b.data)
    if isempty(groups)
        throw(ArgumentError("$b does not contain any non-empty data frames"))
    end
    subbootstrap_parameters = if get(b.parameters, :scheme, :stratified) == :truncated
        (; indices=true, truncate=minimum(length.(groups)) *
            get(b.parameters, :truncationratio, 1))
    else
        (; indices=true)
    end
    train_groups, test_groups =
        mapreduce(traintest,
                  function ((train_acc, test_acc), (train2, test2))
                      push!(train_acc, train2)
                      push!(test_acc, test2)
                      (train_acc, test_acc)
                  end,
                  Bootstrap.(groups, Ref(merge(b.parameters, subbootstrap_parameters))),
                  init = (Vector{Int}[], Vector{Int}[]))
    train_rows, test_rows = reduce(vcat, train_groups), reduce(vcat, test_groups)
    if get(b.parameters, :ensurenonempty, true) && 0 in length.((train_rows, test_rows))
        # Re-run if asked for
        traintest(b)
    elseif get(b.parameters, :indices, false)
        train_rows, test_rows
    else
        parent = first(b.data).parent
        view(parent, train_rows, :), view(parent, test_rows, :)
    end
end

# ------------
# K-fold
# ------------

mutable struct KFold{D} <: AbstractSampler{D}
    const data::D
    const k::Int
    const folds::Vector
    currentfold::Int
    const parameters::NamedTuple
    function KFold(data::D, k::Int, parameters::NamedTuple) where {D}
        folds = kfold(data, k, parameters)
        new{D}(data, k, folds, 0, parameters)
    end
end

KFold(data, k::Int; parameters...) =
    KFold(data, k, NamedTuple(parameters))

function traintest(k::KFold)
    k.currentfold += 1
    if k.currentfold > length(k.folds)
        if get(k.parameters, :cycle, false)
            k.currentfold = 1
        else
            throw(ArgumentError("$k has produced all unique folds, and cycle is false."))
        end
    end
    train = reduce(vcat, k.folds[setdiff(axes(k.folds, 1), k.currentfold)])
    test = k.folds[k.currentfold]
    if get(k.parameters, :invert, false)
        test, train
    else
        train, test
    end
end

function kfold(r::OrdinalRange, k::Integer, ::NamedTuple)
    sr = shuffle(r)
    kindicies =  round.(Int, range(1, length(r)+1, k+1)) |>
        cuts -> UnitRange.(cuts[1:end-1], cuts[2:end] .- 1)
    getindex.(Ref(sr), kindicies)
end

function kfold(data::Union{<:AbstractDataFrame, <:AbstractMatrix}, k::Integer, ::NamedTuple)
    getindex.(Ref(data), kfold(1:size(data, 1), k), :)
end

function kfold(gdf::GroupedDataFrame, k::Integer, params::NamedTuple)
    kgroups = kfold(1:length(gdf), k)
    group_rows = zip(gdf.starts, gdf.ends) .|> r -> r[1]:r[2] |>
        if get(params, :singular, false) rand else identity end
    krows = getindex.(Ref(group_rows), kgroups) .|>
        if get(params, :singular, false) identity else collect ∘ Iterators.flatten end
    if get(params, :indices, false)
        krows
    else
        getindex.(Ref(gdf.parent), krows, :)
    end
end

function kfold(gdfs::Vector{<:GroupedDataFrame}, k::Integer, params::NamedTuple)
    if !all(gdf -> gdf.parent ===(first(gdfs).parent), gdfs)
        throw(ArgumentError("GroupedDataFrames do not all share the same parent"))
    end
    groups = filter(g -> length(g) > 0, gdfs)
    if isempty(groups)
        throw(ArgumentError("No non-empty data frames provided"))
    end

    krows = [kfold(g, k, merge(params, (; indices=true)))
             for g in groups]
    if get(params, :scheme, :stratified) == :truncated
        for i in 1:k
            maxlen = mapreduce(length, min, getindex.(krows, k)) *
                get(params, :truncationratio, 1)
            for krow in krows
                krow[k] = krow[k][1:min(length(krow[k], maxlen))]
            end
        end
    end

    krows = reduce(vcat, krows)
    if get(params, :indices, false)
        krows
    else
        view.(Ref(first(groups).parent), krows, :)
    end
end
