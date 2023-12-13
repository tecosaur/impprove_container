using AverageShiftedHistograms

using Makie

import IMPPROVE.comparemodels

function comparemodels(category1::String, category2::String; modeldir::AbstractString=IMPPROVE.modeldir(),
                       kwargs...)
    list1 = IMPPROVE.listmach(; category=category1, modeldir)
    list2 = IMPPROVE.listmach(; category=category2, modeldir)
    listcommon = list1 ∩ list2
    comparemodels(listcommon,
                  getindex.(IMPPROVE.getmach.(listcommon, :scores; category=category1, modeldir), :scores),
                  getindex.(IMPPROVE.getmach.(listcommon, :scores; category=category2, modeldir), :scores))
end

function comparemodels(labels::Vector{String}, scores1::Vector{Dict{Symbol, Any}}, scores2::Vector{Dict{Symbol, Any}};
                       score::Symbol=:πα⁻¹5000)
    scatter(getindex.(scores1, score), getindex.(scores2, score))
end
