using IMPPROVE, MLJ, MLJDecisionTreeInterface, JLSO, Scratch
using CodecZlib, VariantCallFormat, BioSequences
using ProgressMeter
using DataFrames, CSV

include("apply_lib.jl")
include("render_lib.jl")

haskey(ENV, "_JL_DOCKER_PRERUN") && exit(0)

const SETUP =
    (hpos = Ref{Union{Nothing, Vector{Int}}}(nothing),
     modelsets = Ref(String[]))

if (mflag = readflag!(ARGS, ("-m", "--models"))) |> !isnothing
    availible = filter(isdir, readdir("/models", join=true))
    asked = String.(split(mflag, ','))
    selected = String[]
    for model in asked
        if any(model .== basename.(availible))
            push!(selected, model)
        elseif sum(occursin.(model, availible)) == 1
            push!(selected, basename(availible[findfirst(occursin.(model, availible))]))
        elseif sum(occursin.(lowercase(model), lowercase.(availible))) == 1
            push!(selected, basename(availible[findfirst(occursin.(lowercase(model), lowercase.(availible)))]))
        elseif sum(endswith.(availible, model)) == 1
            push!(selected, basename(availible[findfirst(endswith.(availible, model))]))
        elseif sum(endswith.(lowercase.(availible), lowercase(model))) == 1
            push!(selected, basename(availible[findfirst(endswith.(lowercase.(availible), lowercase(model)))]))
        else
            @warn "No model called '$model' could be identified from the $(length(availible)) availible\
                   \n -$(join(basename.(availible), "\n -"))"
        end
    end
    SETUP.modelsets[] = selected
else
    SETUP.modelsets[] = map(basename, filter(isdir, readdir("/models", join=true)))
end

if !isempty(ARGS)
    SETUP.hpos[] = parsehpo.(ARGS)
end

IMPPROVE.set_modeldir!("/models")

VCF_FILES, HPO_MAP = let files = readdir("/vcfs", join=true)
    vcfs = filter(f -> endswith(f, ".vcf") || endswith(f, ".vcf.gz"),
           files)
    vcfs,
    if "/vcfs/hpo.map" in files
        open("/vcfs/hpo.map") do io
            mapping = Dict{String, Vector{Int}}()
            for line in eachline(io)
                vcf, all_terms = split(line, limit=2)
                if string("/vcfs/", vcf) ∉ vcfs
                    @warn "VCF in hpo.map not present ($vcf), skipping"
                    continue
                end
                terms = map(parsehpo, split(all_terms))
                mapping[vcf] = terms
            end
            mapping
        end
    else
        Dict{String, Vector{Int}}()
    end
end

if isempty(HPO_MAP) && isnothing(SETUP.hpos[])
    @error "No HPO terms specified"
    exit(1)
end

if isempty(VCF_FILES)
    @error "No VCF files found in /vcfs, aborting!"
    exit(1)
end

@info "Loading $(length(VCF_FILES)) VCFs"

VCF_TABLES = DataFrame[]
VCF_ANNOT = Vector{DataFrame}[]

for vcf in VCF_FILES
    table = IMPPROVE.vcfdata(vcf)
    if !isempty(table)
        push!(VCF_TABLES, table)
    else
        @warn "$vcf is empty, skipping"
        deleteat!(VCF_FILES, length(VCF_TABLES)+1)
    end
end

if isempty(VCF_FILES)
    error("No variants to process, exiting.")
end

@info "Annotating the VCFs"

for i in eachindex(VCF_TABLES)
    @debug "Expression sources: $(join(map(guess_base_expr, SETUP.modelsets[]), ", "))"
    push!(VCF_ANNOT, [
        vcfannotate(VCF_TABLES[i], base, expr)
        for (base, expr) in map(guess_base_expr, SETUP.modelsets[])])
end

const MACH_SCORE = :πα⁻¹5000
const MIN_GOOD_SCORE = 0.4
const HEATMAP_TOPN = 50
const HGNC = d"HGNC"

for (i, vcf) in enumerate(VCF_FILES)
    @info "Running predictions on $(basename(vcf))"
    mkpath(joinpath("/predictions", basename(vcf)))
    isnothing(VCF_TABLES[i]) && continue
    annotdata = VCF_ANNOT[i]
    hpos = if haskey(HPO_MAP, basename(vcf))
        HPO_MAP[basename(vcf)]
    elseif !isempty(SETUP.hpos[])
        SETUP.hpos[]
    else
        @warn "No HPO terms attached to $vcf, skipping"
        continue
    end
    unionrows = Dict{UInt64, NamedTuple}()
    for annot in annotdata
        for row in Tables.namedtupleiterator(annot)
            rowid = hash((row.chromosome, row.location, row.change...))
            get!(function ()
                     (; ensembl_id, chromosome, location) = row
                     ref, alt = row.change
                     gene_idx = findfirst(
                         g -> !ismissing(g) && g == ensembl_id,
                         HGNC.ensembl_gene_id)
                     gene_name, gene_symbol = if !isnothing(gene_idx)
                         HGNC[gene_idx, [:name, :symbol]]
                     else
                         missing, missing
                     end
                     (; ensembl_id, gene_symbol, gene_name, chromosome, location, ref, alt)
                 end,
                 unionrows, rowid)
        end
    end
    resdf = DataFrame(collect(values(unionrows)))
    sort!(resdf, [:chromosome, :location], by=[c -> chrom_indexmap[c], identity])
    resdf_rowids = map(function ((; chromosome, location, ref, alt),)
                           change = ref => alt
                           hash((; chromosome, location, change))
                       end,
                       Tables.namedtupleiterator(resdf))
    all_preds = fill(NaN, size(resdf, 1), length(hpos) * length(SETUP.modelsets[]))
    mach_scores = Float64[]
    importances = NamedTuple[]
    labels = String[]
    for (m, (annotdf, model)) in enumerate(zip(annotdata, SETUP.modelsets[]))
        X = select(annotdf, Not([:chromosome, :location, :change, :ensembl_id]))
        df_rowids = map(((; chromosome, location, change),) -> hash((; chromosome, location, change)),
                        Tables.namedtupleiterator(annotdf))
        res_rows = filter(!isnothing,
                          map(id -> findfirst(==(id), resdf_rowids),
                              df_rowids))
        for (i, hpo) in enumerate(hpos)
            @debug "Predicting $hpo with $model X$(size(X))"
            mach = IMPPROVE.getmach(string(hpo), :machine, :scores, category = model)
            push!(mach_scores, mach[:scores][MACH_SCORE])
            preds = if isnothing(mach[:machine]) # CADD/AM
                push!(importances, (; ))
                X[!, end] # The last column should be the CADD/AM score
            else
                push!(importances, mach[:machine].fitresult.importances |> NamedTuple)
                pdf.(MLJ.predict(mach[:machine], X), true)
            end
            push!(labels, "$model " * "HP:" * lpad(hpo, 7, '0'))
            resdf[!, labels[end]] = zeros(Float64, size(resdf, 1))
            all_preds[res_rows, (m-1) * length(hpos) + i] = resdf[res_rows, labels[end]] = preds
        end
    end

    CSV.write(joinpath("/predictions", basename(vcf), "all-predictions.csv"),
              resdf)

    CSV.write(joinpath("/predictions", basename(vcf), "model-scores.csv"),
              DataFrame("model" => labels,
                        "$MACH_SCORE score" => mach_scores))

    for (set, imps, labs) in zip(SETUP.modelsets[],
                                 Iterators.partition(importances, length(hpos)),
                                 Iterators.partition(labels, length(hpos)))
        isempty(first(imps)) && continue
        impdf = DataFrame(imps)
        impdf.model = labs
        impdf′ = unstack(stack(impdf), :variable, :model, :value)
        sort!(impdf′, :variable)
        CSV.write(joinpath("/predictions", basename(vcf), "$set-variable-importances.csv"),
                impdf′)
    end

    CSV.write(joinpath("/predictions", basename(vcf), "hpo-descriptions.csv"),
              select(filter(:hpo => h -> h ∈ hpos, d"HPO descriptions"),
                     :hpo => ByRow(h -> "HP:" * lpad(h, 7, '0')) => :hpo,
                     :description))

    hpo_names = let hdesc = d"HPO descriptions"
        [let idx = findfirst(==(hpo), hdesc.hpo)
             model_shortname(model) * ": " * if !isnothing(idx)
                 hdesc.description[idx]
             else
                 "HP:" * lpad(hpo, 7, '0')
             end
         end
         for model in SETUP.modelsets[] for hpo in hpos]
    end

    variant_labels =
        string.(resdf.gene_symbol, " (chr", resdf.chromosome, ':', resdf.location, ' ', resdf.ref, '>', resdf.alt, ')')
    wsum_order = sortperm(replace(all_preds, NaN => 0.0) * mach_scores, rev=true)
    topn = min(size(resdf, 1), HEATMAP_TOPN)
    save(joinpath("/predictions", basename(vcf), "top-$topn-preds.png"),
         pred_headmap(all_preds[wsum_order[1:topn], :], clust=false,
                      title = "Heatmap of top $topn predictions for $(basename(vcf))",
                      xlabel = "Variant",
                      xticks = (1:topn, variant_labels[wsum_order[1:topn]]),
                      ylabel = "Model",
                      yticks = (1:length(hpo_names), hpo_names),
                      xticklabelrotation = 1,
                      ),
         pt_per_unit = 4.0)
    save(joinpath("/predictions", basename(vcf), "all-preds.png"),
         pred_headmap(all_preds[wsum_order, :], clust=false,
                      title = "Heatmap of predictions for $(basename(vcf))",
                      xlabel = "Variant",
                      xticks = (1:size(all_preds, 1), variant_labels[wsum_order]),
                      ylabel = "Model",
                      yticks = (1:length(hpo_names), hpo_names),
                      xticklabelrotation = 1,
                      ),
         pt_per_unit = 3.0)

    # Good mean
    good_selection = mach_scores .>= MIN_GOOD_SCORE
    good_preds = all_preds[:, good_selection] * mach_scores[good_selection]
    @info "$(sum(mach_scores .>= MIN_GOOD_SCORE)) / $(length(mach_scores)) models have a $MACH_SCORE score ≥ $MIN_GOOD_SCORE (and are considered 'good')"
    good_preds_df = select(resdf, :ensembl_id, :gene_name, :chromosome, :location, :ref, :alt)
    good_preds_df.prediction = good_preds
    CSV.write(joinpath("/predictions", basename(vcf), "wmean-good-prediction.csv"),
              good_preds_df)
end

@info "Done"
