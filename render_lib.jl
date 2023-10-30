using CairoMakie, Distances, Clustering

function pred_headmap(predmat; clust::Bool=true, axargs...)
    if clust
        x_distmat = pairwise(Euclidean(), predmat')
        y_distmat = pairwise(Euclidean(), predmat)
        xclus = hclust(x_distmat).order
        yclus = hclust(y_distmat).order
    else
        xclus = axes(predmat, 1)
        yclus = axes(predmat, 2)
    end
    fig = Figure(resolution = ((20, 30) .* size(predmat) .+ (450, 400)))
    ax = Axis(fig[1,1]; axargs...)
    hmp = heatmap!(ax, predmat[xclus, yclus])
    Colorbar(fig[1,2], hmp, label = "Prediction score")
    fig
end
