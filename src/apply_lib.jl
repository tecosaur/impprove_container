using GeneVariantsMinimalData, DataToolkit
using Base.Threads
using StatsBase

function readflag!(args, forms)
    index = findfirst(a -> a ∈ forms, args)
    if !isnothing(index)
        val = if index < length(args) args[index+1] end
        deleteat!(args, ifelse(isnothing(val), [index], [index, index+1]))
        return val
    end
    index = findfirst(a -> any(f -> startswith(a, f * '='), forms), args)
    if !isnothing(index)
        val = last(split(args[index], '=', limit=2))
        deleteat!(args, index)
        return val
    end
    index = findfirst(a -> any(f -> length(f) == 2 && startswith(a, f), forms), args)
    if !isnothing(index)
        val = args[index][3:end]
        deleteat!(args, index)
        return val
    end
end

const hpodescs = haskey(ENV, "_JL_DOCKER_PRERUN") ||
    d"HPO descriptions"

function parsehpo(h::AbstractString)
    if all(c -> c ∈ '0':'9', strip(h))
        parse(Int, h)
    elseif startswith(lowercase(h), "hp") && h[end] ∈ '0':'9'
        parse(Int, replace(lowercase(h), r"hp:?" => ""))
    else
        regularise(str) = replace(lowercase(str), r"[^a-z0-9]" => "")
        descs = regularise.(hpodescs.description)
        hr = regularise(h)
        inds = findall(d -> startswith(d, hr), descs)
        length(inds) == 1 && return hpodescs.hpo[inds[1]]
        inds = findall(d -> occursin(hr, d), descs)
        length(inds) == 1 && return hpodescs.hpo[inds[1]]
        hparts = regularise.(split(h))
        function loosematch(parts, whole)
            start, i = ncodeunits(whole), 1
            for p in parts
                r = findnext(p, whole, i)
                isnothing(r) && return
                i = first(r)
                start = min(start, i)
            end
            start
        end
        lscores = filter(!isnothing ∘ last,
                         map(d -> loosematch(hparts, d), descs) |>
                             enumerate |> collect)
        isempty(lscores) && error("Could not identify HPO term for '$h'")
        hpodescs.hpo[first(argmin(last, lscores))]
    end
end

function guess_base_expr(name::String)
    name = lowercase(name)
    base = if !isnothing(match(r"(?:_|\b|^)cadd(?:_|\b|$)", name))
        :CADD
    elseif !isnothing(match(r"(?:_|\b|^)am(?:_|\b|$)", name))
        :AM
    else
        error("Unable to guess the base score for '$name'")
    end
    expr = if endswith(name, "_cadd") || endswith(name, "_am")
        ""
    elseif endswith(name, "_gtex")
        "GTEx v8 gene tissue summary"
    elseif endswith(name, "_tabsap_tissue")
        "Tabula Sapiens Tissue summary"
    elseif endswith(name, "_tabsap_tissue_mean")
        "Tabula Sapiens Tissue mean"
    elseif endswith(name, "_tabsap_tissue_dispersion")
        "Tabula Sapiens Tissue dispersion"
    elseif endswith(name, "_tabsap_celltype")
        "Tabula Sapiens Celltype summary"
    elseif endswith(name, "_tabsap_celltype_mean")
        "Tabula Sapiens Celltype mean"
    elseif endswith(name, "_tabsap_celltype_dispersion")
        "Tabula Sapiens Celltype dispersion"
    else
        error("Unable to guess the expression data for '$name'")
    end
    base, expr
end

const chrom_indexmap = indexmap(vcat(string.(1:22), "X", "Y", "M", "MT"))

function vcfannotate(variants::DataFrame, base::Symbol, expr_d::String)
    @debug "Annotating with" base expr_d
    variants = copy(variants) # protect against mutation
    exons = d"Gencode exons"
    expr = if !isempty(expr_d) read(dataset(expr_d)) end
    rows = Tables.rowtable(variants)
    variants.ensembl_id = missings(String, nrow(variants))
    chroms = Dict{String, UnitRange{Int}}()
    extra_lock = SpinLock()
    extras = NamedTuple{(:chromosome, :location, :change, :ensembl_id),
                        Tuple{String, Int, Pair{Char, Char}, Union{String, Missing}}}[]
    nonexon_removed = Threads.Atomic{Int}(0)
    @threads for i in axes(variants, 1)
        (; chromosome, location, change) = rows[i]
        ensembl_ids = GeneVariantsMinimalData.find_surrounding_exons_genes(location, chromosome, exons, protprefer=true)
        if isempty(ensembl_ids)
            nonexon_removed[] += 1
        else
            variants.ensembl_id[i] = first(ensembl_ids)
        end
        if length(ensembl_ids) > 1
            @lock extra_lock for ensembl_id in ensembl_ids[2:end]
                push!(extras, (; chromosome, location, change, ensembl_id))
            end
        end
    end
    # TODO handle extras
    ndropped = sum(ismissing, variants.ensembl_id)
    if ndropped > 0
        @warn "Dropping $ndropped variants ($(round(100*ndropped/size(variants,1), digits=1))%) \
               as they could not be linked to a gene thorugh the exome."
    end
    dropmissing!(variants, :ensembl_id)
    num_have_ensembl = size(variants, 2)
    if base == :CADD
        getcadd = d"getCADD"
        annotate_cadd!(variants, getcadd)
    elseif base == :AM
        alphamissense = d"AlphaMissense Hg38"
        annotate_alphamissense!(variants, alphamissense)
    else
        error("Invalid base score: $base")
    end
    num_have_score = size(variants, 2)
    if num_have_ensembl > num_have_score
        ndropped = num_have_ensembl - num_have_score
        @warn "Dropping $ndropped variants ($(round(100*ndropped/num_have_ensembl, digits=1))%) \
               as they do not have a $base score."
    end
    if !isnothing(expr)
        leftjoin!(variants, expr, on=:ensembl_id)
    end
    dropmissing!(variants)
    num_have_expr = size(variants, 2)
    if num_have_score > num_have_expr
        ndropped = num_have_score - num_have_expr
        @warn "Dropping $ndropped variants ($(round(100*ndropped/num_have_score, digits=1))%) \
               as they do not have a $base score."
    end
    variants
end

function annotate_cadd!(annot::DataFrame, getcadd::Function)
    annot.cadd = zeros(Float64, size(annot, 1))
    for (i, (; chromosome, location, change)) in
        enumerate(Tables.namedtupleiterator(annot))
        annot.cadd[i] = if first(chromosome) != 'M'
            Base.invokelatest(getcadd, chromosome, location, change)
        else
            NaN
        end
    end
    filter!(:cadd => !isnan, annot)
end

model_shortname(name) =
    replace(name,
        "clinvar_gnomad_cadd" => "CADD",
        "clinvar_gnomad_am" => "AM",
        "_tabsap_tissue" => " + TabSap tissue",
        "_tabsap_celltype" => " + TabSap celltype",
        "_gtex" => " + GTEx",
        "_mean" => " mean",
        "_var" => " variance",
        "_dispersion" => " dispersion")

function annotate_alphamissense!(annot::DataFrame, alphamissense::DataFrame)
    alphamissense_scores = zeros(Float64, size(annot, 1))
    current_chr = ""
    current_chr_range = 0:0
    prog = Progress(size(annot, 1), desc="Annotating: ", color=:light_black)
    for (row, (; chromosome, location, change)) in enumerate(Tables.namedtupleiterator(annot))
        next!(prog)
        if chromosome != current_chr
            current_chr = chromosome
            current_chr_range = findfirst(==(chromosome), alphamissense.chromosome):findlast(==(chromosome), alphamissense.chromosome)
        end
        locmatch = searchsorted(view(alphamissense.location, current_chr_range), location)
        isempty(locmatch) && continue
        snvs = view(view(alphamissense.change, current_chr_range), locmatch)
        changematch = findfirst(snvs .== change)
        isnothing(changematch) && continue
        alphamissense_scores[row] = alphamissense.pathogenicity[
            first(current_chr_range) + first(locmatch) + changematch - 2]
    end
    annot.am_score = alphamissense_scores
    filter!(:am_score => >(0), annot)
end
