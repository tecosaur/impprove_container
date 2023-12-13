module SaveExt

using IMPPROVE

using JLSO
using Scratch
using MLJ
using Dates

import IMPPROVE: set_modeldir!, modeldir, listmach,
    getmach, getmodel, savemach, savemodel

MODEL_DIR = get(ENV, "IMPPROVE_MODEL_DIR", @get_scratch!("models"))

set_modeldir!(path::AbstractString) = global MODEL_DIR = path
modeldir() = MODEL_DIR

function listmach(; category::AbstractString="misc", modeldir::AbstractString=MODEL_DIR)
    replace.(filter(p -> endswith(p, ".jlso"), readdir(joinpath(modeldir, category))),
             ".jlso" => "")
end

function getmach(name::AbstractString, attrs::Symbol...; category::AbstractString="misc",
                 modeldir::AbstractString=MODEL_DIR)
    mpath = joinpath(modeldir, category, name * ".jlso")
    if isfile(mpath)
        info = JLSO.load(mpath, attrs...)
        if (isempty(attrs) || :machine in attrs) && !isnothing(get(info, :machine, nothing))
            restore!(info[:machine])
        end
        info
    else
        error("Saved machine $name does not exist [category=$category, modeldir=$modeldir]")
    end
end

getmodel(name::AbstractString; kwargs...) =
    getmach(name, :machine; kwargs...)[:machine].fitresult

function savemach(mach::Machine, label::Any; params::NamedTuple=(;),
                  scores::Union{Nothing, Dict}=nothing,
                  modeldir::AbstractString=MODEL_DIR)
    mpath = joinpath(modeldir, get(params, :name, "misc"), string(label) * ".jlso")
    JLSO.save(mpath,
              :machine => serializable(mach), :label => label,
              :parameters => params, :scores => scores, :created => Dates.now())
end

function savemodel(model, X, y, label; params::NamedTuple=(;), modeldir::AbstractString=MODEL_DIR)
    mach = machine(model, X, y)
    fit!(mach)
    score = if y isa MLJ.CategoricalArrays.CategoricalVector{Bool} &&
               hasfield(mach.fitresult, :oob)
        scoreproba(mach.fitresult.oob)
    end
    savemach(mach, label; params, score, modeldir)
end

end
