using InlineStrings

import IMPPROVE.visualisemodel

function visualisemodel(data::Dict{Symbol, Any}, hpostats::DataFrame;
                        resolution::Tuple{Int,Int} = (1600, 1000),
                        fixaspectratio::Bool=false, getpopulation::Function)
    fig = Figure(; resolution)

    terms = map(keys(data[:terms]) |> collect) do term
        if any(Base.Fix1(isa, term), (Integer, String))
            (term,)
        elseif term isa Tuple
            term
        elseif term isa Vector
            Tuple(term)
        end end |> sort

    logp = Observable(false)
    log1(x) = log(1 + (ℯ-1)*max(0, x))

    centrep = Observable(false)
    πα⁻¹options = Int[50_000, 20_000, 10_000, 5000, 2000, 1000, 500, 200, 100]
    πα⁻¹thresh = Observable(5000)

    importancesort = Observable(false)

    hpo = Observable{Tuple{Vararg{Int}}}(tuple(first(data[:terms])...))
    hdat_cache = Dict{Tuple{Tuple{Vararg{Int}}, Int, Bool}, NamedTuple}()

    hpatho_cache = Dict{Tuple{Vararg{Int}}, DataFrame}()

    function gen_hdat(hpo::Tuple{Vararg{Int}}, αthresh::Number, importancesort::Bool)
        if !haskey(hpatho_cache, hpo)
            hpatho_cache[hpo] = getpopulation(collect(hpo); data[:parameters]...)
        end
        hpatho = hpatho_cache[hpo]
        @info "" names(hpatho)
        pathogenes = unique(hpatho.hgnc_id)
        preds, stat_importances = if eltype(data["predictions"]).parameters[1] <: Tuple
            data["predictions"][hpo],
            data["statistics"][hpo][:importances]
        else
            data["predictions"][if length(hpo) == 1; hpo[1] else hpo end],
            data["statistics"][if length(hpo) == 1; hpo[1] else hpo end][:importances]
        end
        pathopreds = filter(:gene => g -> g ∈ pathogenes, preds)
        benignpreds = filter(:gene => g -> g ∉ pathogenes, preds)
        rates = VARPPI.scorerates(preds.pathogenic, preds.score)

        if importancesort
            hpatho = select(hpatho, first.(sort(collect(stat_importances), by=last, rev=true)), :)
        end

        threshold = rates[:thresholds][last(VARPPI.πinterpolated(rates[:π], rates[:α], αthresh))]
        pathopreds.predicted = pathopreds.score .> threshold
        pathopreds.correct = pathopreds.pathogenic .== pathopreds.predicted

        geneinfo = combine(groupby(pathopreds, :gene),
                           :correct => mean,
                           :correct => sum,
                           nrow)

        expressioncols = setdiff(names(hpatho), ["pathogenic", "hgnc_id", "cadd"])
        geneexpression = combine(groupby(hpatho, :hgnc_id),
                                expressioncols .=> first .=> expressioncols)
        exprtypes = split.(expressioncols, '_') .|> last |> unique
        sort!(geneexpression, :hgnc_id, by = h -> findfirst(==(h), geneinfo.gene))
        select!(geneexpression, :hgnc_id => :gene, Not(:hgnc_id))
        geneexpression = geneexpression[!, sortperm(names(geneexpression),
                                                    by=c -> findfirst(==(last(split(c, '_'))), ["gene"; exprtypes]))]
        exprimportances = fill(0.0, size(geneexpression, 2) - 1)
        for (var, imp) in stat_importances
            var = replace(var, "_ndispersion" => "dispersion")
            pos = findfirst(==(var), names(geneexpression)[2:end])
            isnothing(pos) || (exprimportances[pos] = imp)
        end

        revfirstnonzero = findfirst(>(0.0), geneinfo.correct_mean |> reverse)
        predline = 1 + size(geneinfo, 1) - something(revfirstnonzero, 1)

        pred_subset = geneexpression[1:predline, :]
        nopred_subset = geneexpression[predline+1:end, :]

        predsplitjs = if predline == 0 || predline == size(geneexpression, 1)
            fill(0.0, size(geneexpression, 2) - 1)
        else
            map(names(geneexpression)[2:end]) do var
                VARPPI.jensen_shannon_split(
                    collect(axes(geneexpression, 1)) .<= predline,
                    [geneexpression[!, var]])
            end
        end

        (; hpatho, geneexpression, geneinfo, exprimportances, exprtypes,
         stat_importances, predline, predsplitjs, preds, pathopreds, benignpreds,
         rates, threshold)
    end

    hpo_debounced = async_latest(hpo, 1)

    hdat = map(hpo_debounced, πα⁻¹thresh, importancesort) do hpo, α⁻¹, imps
        get(hdat_cache, (hpo, α⁻¹, imps),
            let hd = gen_hdat(hpo, 1/α⁻¹, imps)
                hdat_cache[(hpo, α⁻¹, imps)] = hd
                hd
            end)
    end

    centre(m::AbstractMatrix; dims::Int=1) = (m .- mean(m; dims)) ./ std(m; dims)

    expr_value_matrix = @lift(Matrix($hdat.geneexpression[!, Not(:gene)]) .|>
        if $logp log1 else identity end |>
        if $centrep m -> clamp.(centre(m), -50, 50) else identity end |> transpose)

    # Inspector labelling

    function plot_label(plot, idx, _pos)
        hgnc = get_artefact(:hgnc)
        var_idx, gene_idx = idx
        val = plot.color[][idx...]
        # val = if plot != hmap_cadd,
        #     expr_value_matrix[][var_idx, gene_idx]
        # else plot[1][][idx] end
        hgnc_id = hdat[][:geneexpression][gene_idx, 1]
        var_name = names(hdat[][:geneexpression])[var_idx+1]
        string(hgnc_id, " (", [hgnc[hgnc.hgnc_id .== hgnc_id, :name]; ""][1], ")\n",
               sum(hdat[][:hpatho].hgnc_id .== hgnc_id), " pathogenic gene variants\n\n",
               if plot != hmap_cadd
                   var_name * "\n\n"
               else "" end,
               "Value: ", round(val, sigdigits = 4))
    end

    function hvar_importance_label(plot, (idx, _), _pos)
        variable = names(hdat[].geneexpression)[idx+1]
        imp = get(hdat[].stat_importances, variable, 0.0)
        string(variable, "\nImportance: ",
               round(imp, sigdigits=3))
    end

    function hvar_correctness_label(plot, (_, idx), _pos)
        string("Correctly predicted: ",
               hdat[].geneinfo.correct_sum[idx],
               '/',
               hdat[].geneinfo.nrow[idx],
               " (",
               round(100 * hdat[].geneinfo.correct_mean[idx],
                     sigdigits = 3),
               "%)")
    end

    function hvar_label(plot, idx, _pos)
        string(names(hdat[][:geneexpression])[Int(idx)+1],
               ": ",
               round(hdat[].predsplitjs[Int(idx)], sigdigits=3))
    end

    # Plotting

    plotmain = fig[1:5,1] = GridLayout()

    ax_cadd = Axis(plotmain[2,1], width=12, yreversed=true,
                   xlabel="CADD", tellwidth=false, xpanlock=true, xrectzoom=false)
    ax_expression = Axis(plotmain[2,2], yreversed=true,
                         xlabel=@lift("Expression " * join($hdat.exprtypes, ", ", ", and ")))
    ax_importances = Axis(plotmain[1,2], height=12, xaxisposition=:top,
                          xlabel="Variable importance")
    ax_cadd_importance = Axis(plotmain[1,1], height=12, width=12,
                              xpanlock=true, xrectzoom=false,
                              ypanlock=true, yrectzoom=false)
    ax_correctness = Axis(plotmain[2,3], width=12, yreversed=true,
                          xlabel="Correctness", tellwidth=false,
                          xpanlock=true, xrectzoom=false)
    ax_jsdiv = Axis(plotmain[3,2], yreversed=true, height=30,
                    xlabel="JS divergence", yrectzoom=false, ypanlock=true)
    ax_predfreq = Axis(fig[3, 2], xlabel="Variant rank",
                       yzoomlock=true, yrectzoom=false, ypanlock=true)
    ax_corrsum = Axis(fig[5, 2], xlabel="Mean gene correctness")

    colsize!(plotmain, 1, Fixed(12))
    if fixaspectratio
        colsize!(plotmain, 2, Aspect(2, 1))
    else
        on(hdat) do hdat
            colsize!(plotmain, 2,
                    Aspect(2, clamp((size(hdat.geneexpression, 2)-1) /
                                        size(hdat.geneexpression, 1),
                                    0.7, 1.2)))
        end
    end
    colsize!(plotmain, 3, Fixed(12))
    rowsize!(plotmain, 1, Fixed(12))

    for ax in (ax_cadd, ax_expression, ax_correctness,
               ax_importances, ax_cadd_importance, ax_jsdiv, ax_corrsum)
        hidespines!(ax)
        ax == ax_cadd || hideydecorations!(ax)
        ax ∈ (ax_expression, ax_corrsum) || hidexdecorations!(ax, label=false)
    end
    linkyaxes!(ax_expression, ax_cadd, ax_correctness)
    linkxaxes!(ax_expression, ax_importances)
    linkxaxes!(ax_expression, ax_jsdiv)

    hmap_cadd =
        heatmap!(ax_cadd, @lift(hcat(combine(groupby($hdat.hpatho, :hgnc_id), :cadd => maximum => :cadd).cadd)'),
                 colormap = :lighttemperaturemap, inspector_label = plot_label)
    hmap_expression =
        heatmap!(ax_expression, expr_value_matrix,
                 colormap=:thermal, inspector_label = plot_label)
    DataInspector(ax_expression)
    hmap_importances =
        heatmap!(ax_importances, @lift((if $logp log1 else identity end).(hcat($hdat.exprimportances))),
                 colormap=:binary, inspector_label = hvar_importance_label)
    hmap_cadd_importance =
        heatmap!(ax_cadd_importance, @lift(hcat(get($hdat.stat_importances, "cadd", 0.0))),
                 colormap=:bilbao, colorrange=@lift(extrema(values($hdat.stat_importances))),
                 inspector_label = (_, _, _) ->
                     string("CADD\nImportance: ",
                            round(get(hdat[].stat_importances, "cadd", 0.0),
                                  sigdigits=3)))
    hmap_correctness =
        heatmap!(ax_correctness, @lift(hcat($hdat.geneinfo.correct_mean)'), colorrange=(0, 1),
                 colormap=:viridis, inspector_label = hvar_correctness_label)
    bp_jsdiv = barplot!(ax_jsdiv, @lift(axes($hdat.predsplitjs, 1)),
                        @lift($hdat.predsplitjs), color=@lift($hdat.predsplitjs),
                        colormap = :blues, inspector_label = hvar_label)
    xlims!(ax_jsdiv, (extrema(axes(hdat[].predsplitjs, 1)) .+ (-0.5, 0.5))...)

    on(hdat, priority=-1) do _
        reset_limits!(ax_expression)
        reset_limits!(ax_jsdiv, xauto=false)
    end

    predline = map(hdat) do hdat
        revfirstnonzero = findfirst(>(0.0), hdat.geneinfo.correct_mean |> reverse)
        1 + size(hdat.geneinfo, 1) - something(revfirstnonzero, 1)
    end

    hlines!.((ax_expression, ax_cadd, ax_correctness),
            map(hdat -> [0.5 + hdat.predline], hdat) |> Ref,
            color = :white, linewidth = 2,
            inspectable = false)

    correct_means = map(hdat) do hdat
        means = hdat.geneinfo.correct_mean
        if all(==(0.0), means)
            means[end] = nextfloat(0.0) # hist has a problem with all equal values
        end
        means
    end
    hist!(ax_corrsum, correct_means,
             bins=15, #color=:transparent,
             # strokecolor=:royalblue, strokewidth=2,
             inspectable = false)
    xlims!(ax_corrsum, 0, 1)
    on(hdat, priority=-1) do _
        reset_limits!(ax_corrsum, xauto=false)
    end

    allpredspos = @lift(
        findall($hdat.preds.pathogenic .== true))

    hist!(ax_predfreq, allpredspos,
          bins=@lift(round(Int, maximum($allpredspos)/10)))
    # density!(ax_predfreq, allpredspos, npoints=1000,
    #          strokewidth = 2, color = :transparent)
    vlines!(ax_predfreq, allpredspos,
            ymin=0, ymax=0.2, color=:orange)
    vlines!(ax_predfreq,
            @lift([sum($hdat.preds.score .>= $hdat.threshold) + 0.5]),
            color=:red)
    hlines!([0], color=:grey)
    xlims!(ax_predfreq, low = 0)
    last_hpo::Tuple{Vararg{Int}} = (-1,)
    on(hdat, update=true) do _
        if last_hpo != hpo[]
            reset_limits!(ax_predfreq, yauto=false)
            yauto = Makie.yautolimits(ax_predfreq)
            ylims!(ax_predfreq, -0.3*yauto[2], yauto[2])
            last_hpo = hpo[]
        end
    end

    colsize!(fig.layout, 1, Relative(0.7))
    rowsize!(fig.layout, 3, Relative(0.16))
    rowsize!(fig.layout, 5, Relative(0.08))

    # Interface

    interface = fig[1,2] = GridLayout()

    local hpoindex = Observable(first(axes(hpos, 1)))

    if length(data["statistics"]) > 1
        hpo_back = Button(interface[4,1], label="◀", width=30)
        hpo_forward = Button(interface[4,2], label="▶", width=30)

        hpo_slider = Slider(interface[4,3], range=axes(hpos, 1))
        connect!(hpoindex, hpo_slider.value)

        hpo_text = Textbox(interface[4,4], placeholder = "HP:?", width=80)
        hpo_text_do = Button(interface[4,5], label="Go")

        on(idx -> hpo[] = hpos[idx], hpoindex)
        on(_ -> hpoindex[] = max(1, hpoindex[] - 1), hpo_back.clicks)
        on(_ -> hpoindex[] = min(length(hpos), hpoindex[] + 1), hpo_forward.clicks)

        on(hpo_text_do.clicks) do _
            hpo_text.stored_string[] = hpo_text.displayed_string[]
        end

        on(hpo_text.stored_string) do stored_string
            hpo_id = something(tryparse(Int, stored_string), hpo[])
            idx = findfirst(==(tuple(hpo_id)), hpos)
            isnothing(idx) || (hpoindex[] = idx)
        end
    end

    settings_grid = interface[5,1:5] = GridLayout()

    log_switch = Toggle(settings_grid[1,1])
    Label(settings_grid[1,2], "log", halign=:left)
    connect!(logp, log_switch.active)

    centre_switch = Toggle(settings_grid[1,3])
    Label(settings_grid[1,4], "centre", halign=:left)
    connect!(centrep, centre_switch.active)

    imps_switch = Toggle(settings_grid[1,5])
    Label(settings_grid[1,6], "imp. sort", halign=:left)
    connect!(importancesort, imps_switch.active)

    Label(interface[2,1:3], @lift(join([string("HP:", lpad(hpoid, 7, '0'))
                                        for hpoid in vcat(collect($hpo))], " ")),
          textsize=24, halign=:left, tellwidth=false)
    if length(hpos) > 1
        Label(interface[2,4:5], @lift(string('[', $hpoindex, '/', length(hpos), ']')),
              textsize=24, halign=:left, tellwidth=false, tellheight=false)
    end
    Label(interface[3,1:5], @lift(hpostats[hpostats.hpo .== $hpo, :desc]),
          halign=:left, tellwidth=false)

    rowsize!(interface, 1, Fixed(0))
    rowgap!(interface, 2, Fixed(2))
    rowgap!(interface, 4, Fixed(8))
    rowgap!(interface, 4, Fixed(12))
    colgap!(fig.layout, 1, Fixed(40))

    # ROC and PRC plots, threshold info and setting

    thresholdinfo = fig[2,2] = GridLayout()

    closestthresholdidx = map(hdat) do hdat
        argmin(abs.(hdat.threshold .- hdat.rates[:thresholds]))
    end

    Label(thresholdinfo[1,1:2], text=
        @lift(string("Threshold: ", round($hdat.threshold, digits=4),
                     ", α: ", round($hdat.rates[:α][$closestthresholdidx], sigdigits=2),
                     ", π: ", round($hdat.rates[:π][$closestthresholdidx], sigdigits=2))),
          tellwidth=false)

    ax_roc = Axis(thresholdinfo[2,1], xlabel="α", ylabel="π",
                  limits=(-0.1, 1.1, -0.1, 1.1), aspect=AxisAspect(1))
    ax_prc = Axis(thresholdinfo[2,2], xlabel="precision", ylabel="recall",
                  limits=(-0.1, 1.1, -0.1, 1.1), aspect=AxisAspect(1))

    function inspect_xylabels(xlab::String, ylab::String; roundargs...)
        function (plot, idx, _)
            x, y = plot[1][][idx]
            string(xlab, ": ", round(x; roundargs...),
                   "\n", ylab, ": ", round(y; roundargs...))
        end
    end

    lines!(ax_roc, @lift(zip($hdat.rates[:α], $hdat.rates[:π]) |> collect),
           inspectable = false)
    scatter!(ax_roc, @lift(zip($hdat.rates[:α], $hdat.rates[:π]) |> collect),
             color=@lift($hdat.rates[:thresholds]),
             inspector_label = inspect_xylabels("α", "π", digits=6))

    lines!(ax_prc, @lift(zip($hdat.rates[:precision], $hdat.rates[:recall]) |> collect),
           inspectable = false)
    scatter!(ax_prc, @lift(zip($hdat.rates[:precision], $hdat.rates[:recall]) |> collect),
             color=@lift($hdat.rates[:thresholds]),
             inspector_label = inspect_xylabels("precision", "recall", digits=6))

    scatter!(ax_roc, @lift([($hdat.rates[:α][$closestthresholdidx],
                             $hdat.rates[:π][$closestthresholdidx])]),
             color=:red, inspector_label = inspect_xylabels(
                 "current α", "current π", digits=6))
    scatter!(ax_prc, @lift([($hdat.rates[:precision][$closestthresholdidx],
                            $hdat.rates[:recall][$closestthresholdidx])]),
             color=:red, inspector_label = inspect_xylabels(
                 "current precision", "current recall", digits=6))

    thresh_set_group = thresholdinfo[3,1:2] = GridLayout()

    Label(thresh_set_group[1,1], text=@lift(string("Max π such that α ≤ 1/",
                                                   $πα⁻¹thresh, ".")),
          tellwidth=false)
    Label(thresh_set_group[1,2], text="max(α⁻¹):", tellwidth=false)
    thresh_slider = Slider(thresh_set_group[1,3], range=axes(πα⁻¹options, 1),
                           startvalue=findfirst(==(πα⁻¹thresh[]), πα⁻¹options),
                           width=120)
    colgap!(thresh_set_group, 2, 10)

    on(async_latest(thresh_slider.value, 1)) do idx
        πα⁻¹thresh[] = πα⁻¹options[idx]
    end

    summary_table = map(hdat) do (; pathopreds, benignpreds, threshold)
        npatho, nbenign = size(pathopreds, 1), size(benignpreds, 1)
        predpatho = sum(pathopreds.score .> threshold)
        predbenign = sum(benignpreds.score .> threshold)
        table = [
            "(Variants)" "Pathogenic" "Benign";
            "Total" npatho nbenign;
            "Correct" predpatho nbenign-predbenign;
            "Incorrect" npatho-predpatho predbenign;
            "Accuracy" #=
            =# string(round(100*predpatho/npatho, digits=1), "%") #=
            =# string(round(100-100*predbenign/nbenign, digits=3), "%")
        ] .|> string
        colwidths = [maximum(length.(string.(vals))) for vals in eachcol(table)]
        stable = map(enumerate(eachrow(table))) do (rownum, row)
            map(enumerate(row)) do (colnum, val)
                " " * lpad(val, colwidths[colnum]) *
                    if colnum == size(table, 2) " " else "" end
            end
        end
        insert!(stable, 2, ['━'^(1 + sum(colwidths) + length(colwidths))])
        join(join.(stable, ""), "\n")
    end

    Label(fig[4,2], text=summary_table,
          font = if Sys.isapple()
              "SF Mono"
          elseif Sys.iswindows()
              "Consolas"
          else
              "DejaVu Sans Mono"
          end,
          tellwidth=false)

    fig
end

function visualisemodel(name::String, hpostats::DataFrame; kwargs...)
    kw_load = NamedTuple{Tuple(filter(k -> k ∉ (:resolution, :fixedpointnumbers, :getpopulation),
                                      keys(kwargs)))}(NamedTuple(kwargs))
    kw_vis = NamedTuple{Tuple(filter(k -> k ∈ (:resolution, :fixedpointnumbers, :getpopulation),
                                     keys(kwargs)))}(NamedTuple(kwargs))
    mach = loadfromcache(name, :label, :parameters, :scores, :machine; kw_load...)
    data = Dict{Symbol, Any}(:parameters => mach.parameters,
                             :terms => [mach.label],
                             :results => [mach.machine.fitresult])
    visualisemodel(data, hpostats; kw_vis...)
end

visualisemodel(data::Dict{Symbol, Any}; kwargs...) =
    visualisemodel(data; kwargs...)
