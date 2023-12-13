using Random

import MLJBase.train_test_pairs

asrng(::Nothing) = Random.GLOBAL_RNG
asrng(n::Int) = Xoshiro(n)
asrng(rng::AbstractRNG) = rng
const PotentialRNG = Union{Nothing, Int, <:AbstractRNG}

struct Bootstrap <: ResamplingStrategy
    rng::AbstractRNG
    Bootstrap(rng::PotentialRNG) =
        new(asrng(rng))
end

Bootstrap() = Bootstrap(nothing)

function train_test_pairs(b::Bootstrap, rows::AbstractVector)
    train = rand(b.rng, axes(rows, 1), length(rows))
    test = setdiff(axes(rows, 1), train)
    @views [(rows[train], rows[test])]
end

function train_test_pairs(b::Bootstrap, rows::Union{<:AbstractMatrix, <:AbstractDataFrame})
    train, test = train_test_pairs(b, axes(rows, 1)) |> first
    @views [(rows[train, :], rows[test, :])]
end

struct TruncatedBootstrap <: ResamplingStrategy
    rng::AbstractRNG
    n::Int
    discard::Bool
    TruncatedBootstrap(rng::PotentialRNG, n::Int, discard::Bool=false) =
        new(asrng(rng), n, discard)
end

TruncatedBootstrap(n::Int, discard::Bool=false) =
    TruncatedBootstrap(nothing, n, discard)

function train_test_pairs(tb::TruncatedBootstrap, rows::AbstractVector)
    train = rand(tb.rng, axes(rows, 1), min(tb.n, length(rows)))
    test = setdiff(axes(rows, 1), train)
    if tb.discard && tb.n < length(rows)
        test = setdiff(test, rand(tb.rng, axes(rows, 1), length(rows) - tb.n))
    end
    @views [(rows[train], rows[test])]
end

function train_test_pairs(tb::TruncatedBootstrap, rows::Union{<:AbstractMatrix, <:AbstractDataFrame})
    train, test = train_test_pairs(tb, axes(rows, 1)) |> first
    @views [(rows[train, :], rows[test, :])]
end

struct Stratified{S} <: ResamplingStrategy where { S <: ResamplingStrategy }
    rng::AbstractRNG
    sampler::S
    groups::Vector{<:AbstractVector}
    train::Symbol
    test::Symbol
    Stratified(rng::PotentialRNG, sampler::S, groups::Vector{<:AbstractVector},
               train::Symbol, test::Symbol) where { S <: ResamplingStrategy} =
                   new{S}(asrng(rng), sampler, groups, train, test)
end

Stratified(rng::PotentialRNG, sampler::ResamplingStrategy, groups::Vector{<:AbstractVector}; train::Symbol=:group, test::Symbol=:group) =
    Stratified(rng, sampler, groups, train, test)

Stratified(rng::PotentialRNG, sampler::ResamplingStrategy, gdf::GroupedDataFrame, args...; kwargs...) =
    Stratified(rng, sampler,
               getindex.(Ref(gdf.idx), UnitRange.(gdf.starts, gdf.ends)),
               args...; kwargs...)

Stratified(rng::PotentialRNG, sampler::ResamplingStrategy, df::DataFrame, groupcolumn::Union{<:AbstractString, Symbol},
           args...; kwargs...) =
    Stratified(rng, sampler, groupby(df, groupcolumn), args...; kwargs...)

Stratified(samp::ResamplingStrategy, args...; kwargs...) =
    Stratified(nothing, samp, args...; kwargs...)

function train_test_pairs(s::Stratified, rows::AbstractVector)
    grab_rows(groupinds, mode) =
        if mode == :singular
            rand.(s.rng, s.groups[groupinds])
        elseif mode == :samesingular
            rand.(s.rng, s.groups)[groupinds]
        elseif mode == :group
            Iterators.flatten(s.groups[groupinds] .|> collect) |> collect
        else
            throw(ArgumentError("$s contains invalid group mode: $(mode)"))
        end
    train_test_sets = train_test_pairs(s.sampler, axes(s.groups, 1))
    while any(isempty, first.(train_test_sets))
        train_test_sets = train_test_pairs(s.sampler, axes(s.groups, 1))
    end
    [@views (rows[grab_rows(train_grpinds, s.train)],
             rows[grab_rows(test_grpinds, s.test)])
     for (train_grpinds, test_grpinds) in train_test_sets]
end

train_test_pairs(s::Stratified, rows::Union{<:AbstractMatrix, <:AbstractDataFrame}) =
    [@views (rows[train, :], rows[test, :])
     for (train, test) in train_test_pairs(s, axes(rows, 1))]

struct CompositeSampler{S} <: ResamplingStrategy
    subsamplers::S
    CompositeSampler(samplers::Tuple{Vararg{ResamplingStrategy}}) =
        new{typeof(samplers)}(samplers)
    CompositeSampler(samplers::ResamplingStrategy...) =
        new{typeof(samplers)}(samplers)
end

train_test_pairs(s::CompositeSampler, rows::AbstractVector) =
    map([zip(pairs...) for pairs in
             zip(train_test_pairs.(s.subsamplers, Ref(rows))...)]) do (trains, tests)
        @views (rows[Iterators.flatten(trains) |> collect],
                 rows[Iterators.flatten(tests) |> collect])
    end

train_test_pairs(s::CompositeSampler, rows::Union{<:AbstractMatrix, <:AbstractDataFrame}) =
    [@views (rows[train, :], rows[test, :])
     for (train, test) in train_test_pairs(s, axes(rows, 1))]

function twolevelbootstrap(categories::Vector{SubDataFrame}, subgroup::Symbol;
                           maxratio::Number=Inf, train::Symbol=:group, test::Symbol=:group,
                           rng::PotentialRNG=nothing)
    if !all(sdf -> getfield(sdf, :parent) == getfield(first(categories), :parent), categories)
        throw(ArgumentError("All `categories` must be from the same parent DataFrame."))
    end
    catrowgrps = map(categories) do catgrp
        csgrp = groupby(catgrp, subgroup)
        getindex.(Ref(getfield(catgrp, :rows)),
                  getindex.(Ref(csgrp.idx), UnitRange.(csgrp.starts, csgrp.ends)))
    end
    basesampler = if maximum(length.(catrowgrps)) >
        maxratio * minimum(length.(catrowgrps))
        TruncatedBootstrap(rng, round(Int, maxratio * minimum(length.(catrowgrps))))
    else
        Bootstrap(rng)
    end
    CompositeSampler(Stratified.(rng, Ref(basesampler), catrowgrps; train, test) |> Tuple)
end

twolevelbootstrap(gdf::GroupedDataFrame, args...; kwargs...) =
    twolevelbootstrap(Vector{SubDataFrame}(gdf |> collect), args...; kwargs...)

twolevelbootstrap(df::DataFrame, category::Symbol, args...; kwargs...) =
    twolevelbootstrap(groupby(df, category), args...; kwargs...)
