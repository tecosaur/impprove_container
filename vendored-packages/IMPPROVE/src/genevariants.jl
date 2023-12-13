function genesample(df::DataFrame, args...; kwargs...)
    pathogenes = fill(false, size(df, 1))
    map(groupby(df, :ensembl_id) |> collect) do geneset
        if any(geneset[!, :pathogenic])
            pathogenes[getfield(geneset, :rows)] .= true
        end
    end
    dfg = DataFrame(pathogene = pathogenes, gene = df.ensembl_id)
    twolevelbootstrap(dfg, :pathogene, :gene, args...; kwargs...)
end
