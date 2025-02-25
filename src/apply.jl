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
    if !allunique(SETUP.hpos[])
        dups = ["HP:" * lpad(t, 7, '0')
                for t in unique(SETUP.hpos[])
                    if sum(==(t), SETUP.hpos[]) > 1]
        @warn "Removing duplicates of HPO terms: $(join(dups, ", "))"
        unique!(SETUP.hpos[])
    end
end

IMPPROVE.set_modeldir!("/models")

VCF_FILES, HPO_MAP = let files = readdir("/vcfs", join=true)
    itemp = if startswith(first(DEPOT_PATH), tempdir())
        first(DEPOT_PATH)
    else tempdir() end
    vcfs = map(filter(f -> endswith(f, r"\.g?vcf(?:\.gz)?"), files)) do f
        if endswith(f, r"\.gvcf(?:\.gz)?")
            @info "Converting $(basename(f)) to a non-g(vcf)"
            newvcf = itemp * replace(basename(f), r"\.gvcf(?:\.gz)?$" => "-converted.vcf")
            open(newvcf, "w") do io
                for line in eachline(
                    pipeline(`bcftools convert $f
                            --fasta-ref $(abspath("liftover/data/Hg38p14.fa")) --gvcf2vcf`,
                            `bcftools view -e 'ALT[0] == "<NON_REF>"'`))
                    print(io, replace(line, ",<NON_REF>" => ""), '\n')
                end
            end
            newvcf
        else f end
    end
    vcfs,
    if "/vcfs/hpo.map" in files
        vcfnames = Set(map(basename, vcfs))
        open("/vcfs/hpo.map") do io
            mapping = Dict{String, Vector{Int}}()
            for line in eachline(io)
                vcf, all_terms = split(line, limit=2)
                if vcf ∉ vcfnames
                    @warn "VCF in hpo.map not present ($vcf), skipping"
                    continue
                end
                vcf = replace(vcf, r"\.gvcf(?:\.gz)?$" => "-converted.vcf")
                terms = map(parsehpo, split(all_terms))
                if !allunique(terms)
                    dups = ["HP:" * lpad(t, 7, '0')
                            for t in unique(terms) if sum(==(t), terms) > 1]
                    @warn "Removing duplicates of HPO terms: $(join(dups, ", "))"
                    terms = unique(terms)
                end
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

const MACH_SCORE = :πα⁻¹5000
const MIN_GOOD_SCORE = 0.4
const HEATMAP_TOPN = 50
const HEATMAP_MAXN = 500
const HGNC = d"HGNC"
const HDESC = d"HPO descriptions"

const NUM_PROCESSED = Ref(0)

function process_vcf(vcf::String)
    i = (NUM_PROCESSED[] += 1)
    @info "* [$i/$(length(VCF_FILES))] $(basename(vcf))"
    mkpath(joinpath("/predictions", basename(vcf)))
    @info "Loading VCF"
    vcf_table = IMPPROVE.vcfdata(vcf)
    if isempty(vcf_table)
        @warn "$vcf is empty, skipping"
        return
    end
    @info "Annotating (x$(length(SETUP.modelsets[])))"
    @debug "Expression sources: $(join(map(guess_base_expr, SETUP.modelsets[]), ", "))"
    annotdata = [vcfannotate(vcf_table, base, expr)
                 for (base, expr) in map(guess_base_expr, SETUP.modelsets[])]
    hpos = if haskey(HPO_MAP, basename(vcf))
        HPO_MAP[basename(vcf)]
    elseif !isnothing(SETUP.hpos[]) && !isempty(SETUP.hpos[])
        SETUP.hpos[]
    else
        @warn "No HPO terms attached to $vcf, skipping"
        return
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
                     (; chromosome, location, ref, alt, ensembl_id, gene_symbol, gene_name)
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
    all_pred_ranks = fill((0, 0, 0, 0), size(resdf, 1), length(hpos) * length(SETUP.modelsets[]))
    mach_scores = Float64[]
    base_models = Tuple{Int, Symbol}[]
    labels = String[]
    hpo_labels = String[]
    for (m, (annotdf, model)) in enumerate(zip(annotdata, SETUP.modelsets[]))
        X = select(annotdf, Not([:chromosome, :location, :change, :ensembl_id]))
        df_rowids = map(((; chromosome, location, change),) -> hash((; chromosome, location, change)),
                        Tables.namedtupleiterator(annotdf))
        res_rows = filter(!isnothing,
                          map(id -> findfirst(==(id), resdf_rowids),
                              df_rowids))
        importances = NamedTuple[]
        for hpo in hpos
            try
                @debug "Predicting $hpo with $model X$(size(X))"
                mbase, mexpr = guess_base_expr(model)
                if mexpr == ""
                    mbase ∈ last.(base_models) && continue
                    push!(base_models, (length(mach_scores)+1, mbase))
                end
                mach = IMPPROVE.getmach(string(hpo), :machine, :scores, category = model)
                push!(mach_scores, mach[:scores][MACH_SCORE])
                preds = if isnothing(mach[:machine]) # CADD/AM
                    push!(importances, (; ))
                    X[!, end] # The last column should be the CADD/AM score
                else
                    push!(importances, mach[:machine].fitresult.importances |> NamedTuple)
                    pdf.(MLJ.predict(mach[:machine], X), true)
                end
                if mexpr == ""
                    push!(labels, model)
                    push!(hpo_labels, String(mbase))
                else
                    push!(labels, "$model " * "HP:" * lpad(hpo, 7, '0'))
                    idx = findfirst(==(hpo), HDESC.hpo)
                    push!(hpo_labels, model_shortname(model) * ": " * if !isnothing(idx)
                              HDESC.description[idx]
                          else
                              "HP:" * lpad(hpo, 7, '0')
                          end)
                end
                resdf[!, labels[end]] = zeros(Float64, size(resdf, 1))
                all_preds[res_rows, length(labels)] = resdf[res_rows, labels[end]] = preds
                # Rankings
                if !isnothing(mach[:machine])
                    oob_predtrue = pdf.(mach[:machine].fitresult.oob.predicted, true)
                    oob_ispatho = mach[:machine].fitresult.oob.actual
                    oob_patho_preds = sort(oob_predtrue[oob_ispatho])
                    oob_benign_preds = sort(oob_predtrue[.!oob_ispatho])
                    pred_ranks = map(preds) do score
                        p_index = searchsortedlast(oob_patho_preds, score)
                        b_index = searchsortedlast(oob_benign_preds, score)
                        # (tp, fn, tn, fp)
                        (length(oob_patho_preds) - p_index,
                        p_index, b_index,
                        length(oob_benign_preds) - b_index)
                    end
                    all_pred_ranks[res_rows, length(labels)] = pred_ranks
                end
            catch err
                @error sprint(showerror, err)
            end
        end
        (isempty(importances) || isempty(first(importances))) && continue
        impdf = DataFrame(importances)
        impdf.model = labels[end-length(importances)+1:end]
        impdf′ = unstack(stack(impdf), :variable, :model, :value)
        sort!(impdf′, :variable)
        CSV.write(joinpath("/predictions", basename(vcf), "$(first(split(labels[end])))-variable-importances.csv"),
                  impdf′)
    end

    if isempty(labels)
        @warn "No models could be applied"
        return
    end

    all_preds = all_preds[:, 1:length(labels)]

    CSV.write(joinpath("/predictions", basename(vcf), "all-predictions.csv"),
              resdf)

    let confdf = DataFrame(:chromosome => String[], :location => Int[], :ref => Char[], :alt => Char[], :confusion => String[],
                           map(l -> Symbol(l) => Int[], labels)...)
       for (i, row) in enumerate(Tables.namedtupleiterator(resdf))
           for (c, conftype) in enumerate(("tp", "fn", "tn", "fp"))
               confrow = [row.chromosome, row.location, row.ref, row.alt, conftype]
               for j in eachindex(labels)
                   push!(confrow, all_pred_ranks[i, j][c])
               end
               push!(confdf, confrow)
           end
        end
        CSV.write(joinpath("/predictions", basename(vcf), "confusion-matrices.csv"),
                  confdf)
    end

    CSV.write(joinpath("/predictions", basename(vcf), "model-scores.csv"),
              DataFrame("model" => labels,
                        "$MACH_SCORE score" => mach_scores))

    CSV.write(joinpath("/predictions", basename(vcf), "hpo-descriptions.csv"),
              select(filter(:hpo => h -> h ∈ hpos, d"HPO descriptions"),
                     :hpo => ByRow(h -> "HP:" * lpad(h, 7, '0')) => :hpo,
                     :description))

    variant_labels =
        string.(resdf.gene_symbol, " (chr", resdf.chromosome, ':', resdf.location, ' ', resdf.ref, '>', resdf.alt, ')')
    baserescaled_mach_scores = copy(mach_scores)
    for (i, _) in base_models
        baserescaled_mach_scores[i] = baserescaled_mach_scores[i] ./
            -(reverse(extrema(filter(!isnan, all_preds[:, i])))...)
    end
    all_preds_nonan = replace(all_preds, NaN => 0.0)
    wsum_order = sortperm(all_preds_nonan * baserescaled_mach_scores, rev=true)
    topn = min(size(resdf, 1), HEATMAP_TOPN)
    save(joinpath("/predictions", basename(vcf), "top-$topn-preds.png"),
         pred_headmap(all_preds[wsum_order[1:topn], :], clust=false, base_score_inds = first.(base_models),
                      title = "Heatmap of top $topn predictions for $(basename(vcf))",
                      xlabel = "Variant",
                      xticks = (1:topn, variant_labels[wsum_order[1:topn]]),
                      ylabel = "Model",
                      yticks = (1:length(hpo_labels), hpo_labels),
                      xticklabelrotation = 1),
         pt_per_unit = 4.0)
    alln = min(size(resdf, 1), HEATMAP_MAXN)
    save(joinpath("/predictions", basename(vcf), if alln == size(resdf, 1) "all-preds.png" else "top-$alln-preds.png" end),
         pred_headmap(all_preds[wsum_order[1:alln], :], clust=false, base_score_inds = first.(base_models),
                      title = "Heatmap of predictions for $(basename(vcf))",
                      xlabel = "Variant",
                      xticks = (1:size(all_preds, 1), variant_labels[wsum_order]),
                      ylabel = "Model",
                      yticks = (1:length(hpo_labels), hpo_labels),
                      xticklabelrotation = 1),
         pt_per_unit = 3.0)

    # Good mean
    good_selection = mach_scores .>= MIN_GOOD_SCORE
    good_preds_raw = all_preds_nonan[:, good_selection] * baserescaled_mach_scores[good_selection]
    good_preds_weightsum = map(!isnan, all_preds[:, good_selection]) * mach_scores[good_selection]
    good_preds = good_preds_raw ./ good_preds_weightsum
    @info "$(sum(mach_scores .>= MIN_GOOD_SCORE)) / $(length(mach_scores)) models have a $MACH_SCORE score ≥ $MIN_GOOD_SCORE (and are considered 'good')"
    good_preds_df = select(resdf, :chromosome, :location, :ref, :alt, :ensembl_id, :gene_symbol, :gene_name)
    good_preds_df.prediction = good_preds
    CSV.write(joinpath("/predictions", basename(vcf), "wmean-good-prediction.csv"),
              good_preds_df)

    # Mini-readme
    mini_readme = replace(read("/mini-readme.txt", String),
                            "{%infile%}" => basename(vcf))
    write(joinpath("/predictions", basename(vcf), "README.txt"),
          mini_readme)
end

for vcf in VCF_FILES
    process_vcf(vcf)
end

if startswith(first(DEPOT_PATH), tempdir())
    rm(first(DEPOT_PATH), recursive=true, force=true)
end

@info "Done"
