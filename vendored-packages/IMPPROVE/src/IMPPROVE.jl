module IMPPROVE

using DataFrames
using ProgressMeter
using Transducers

using MLJModelInterface
const MMI = MLJModelInterface
using MLJBase

if !isdefined(Base, :get_extension)
    using Requires
end

export ConsensusModel, scoreproba
include("models.jl")
include("scoring.jl")

export MatrixTable
include("matrixtable.jl")

export Bootstrap, TruncatedBootstrap, Stratified, CompositeSampler,
    twolevelbootstrap
include("resampling.jl")
include("genevariants.jl")

include("extfns.jl")

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require MLJDecisionTreeInterface = "c6f25543-311c-4c74-83dc-3ea6d1015661" include("../ext/DecisionTreeExt.jl")
        @require ShapML = "8523bd24-6f66-414b-9e1e-b02db3916d64" include("../ext/ShapExt.jl")
        @require JLSO = "9da8a3cd-07a3-59c0-a743-3fdc52c30d11" begin
            @require Scratch = "6c6a2e73-6563-6170-7368-637461726353" begin
                include("../ext/SaveExt.jl")
                @require AverageShiftedHistograms = "77b51b56-6f8f-5c3a-9cb4-d71f9594ea6e" begin
                    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
                        @require InlineStrings = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48" begin
                            @require Observables = "510215fc-4207-5dde-b226-833fc4488ee2" begin
                                include("../ext/VisualisationExt/VisualisationExt.jl")
                            end end end end end end
        @require VariantCallFormat = "28eba6e3-a997-4ad9-87c6-d933b8bca6c1" begin
            @require InlineStrings = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48" begin
                @require TranscodingStreams = "3bb67fe8-82b1-5028-8e26-92a6c54297fa" begin
                    @require CodecZlib = "944b1d66-785c-5afd-91f1-9de20f533193" begin
                        include("../ext/VCFExt.jl")
                    end end
                @require JLSO = "9da8a3cd-07a3-59c0-a743-3fdc52c30d11" begin
                    include("../ext/VCFSaveExt.jl")
                end end end
    end
end

end
