using CairoMakie, Distances, Clustering

function pred_headmap(predmat; clust::Bool=true, base_score_inds::Vector{Int} = Int[], axargs...)
    if clust
        y_distmat = pairwise(Euclidean(), predmat)
        yclus = hclust(y_distmat).order
    else
        yclus = axes(predmat, 2) |> collect
    end
    filter!(y -> y ∉ base_score_inds , yclus)
    fig = Figure(resolution = ((30, 20) .* size(predmat) .+ (500, 400)))
    axargs = NamedTuple(axargs)
    if !isempty(base_score_inds)
        base_attrs = (:xlabel, :xticks, :xticklabelrotation) ∩ keys(axargs)
        baxargs = NamedTuple{Tuple(base_attrs)}(axargs)
        baxargs = Base.setindex(baxargs, (1:length(base_score_inds), axargs.yticks[2][base_score_inds]), :yticks)
        axargs = NamedTuple{Tuple(setdiff(keys(axargs), base_attrs))}(axargs)
    end
    if haskey(axargs, :yticks) && axargs.yticks isa Tuple
        axargs = Base.setindex(axargs, (axargs.yticks[1][1:length(yclus)], axargs.yticks[2][yclus]), :yticks)
    end
    ax = Axis(fig[1,1]; axargs...)
    !isempty(base_score_inds) && hidexdecorations!(ax)
    hmp = heatmap!(ax, predmat[:, yclus], colorrange=(0.0, 1.0))
    Colorbar(fig[1,2], hmp, label = "Prediction score")
    if !isempty(base_score_inds)
        bax = Axis(fig[2,1]; baxargs...)
        bhmp = heatmap!(bax, predmat[:, base_score_inds])
        rowsize!(fig.layout, 2, Fixed(15))
    end
    fig
end
