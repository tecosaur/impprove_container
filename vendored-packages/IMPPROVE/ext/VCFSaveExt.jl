module VCFSaveExt

using IMPPROVE

using VariantCallFormat
using InlineStrings
using DataFrames
using Logging

using JLSO

import IMPPROVE: vcfcached, vcflog

struct CompositeLogger <: AbstractLogger
    loggers::Vector{<:AbstractLogger}
end

CompositeLogger(loggers...) = CompositeLogger(collect(loggers))

Logging.handle_message(complog::CompositeLogger, level, message, _module, group, id, file, line; kwargs...) =
    for logger in complog.loggers
        Logging.shouldlog(logger, level, _module, group, id) &&
            Logging.handle_message(logger, level, message, _module, group, id, file, line; kwargs...)
    end

Logging.shouldlog(complog::CompositeLogger, level, _module, group, id) =
    any(Logging.shouldlog(logger, level, _module, group, id) for logger in complog.loggers)

Logging.min_enabled_level(complog::CompositeLogger) =
    minimum(map(Logging.min_enabled_level, complog.loggers))

function vcfcached(vcf::Union{VariantCallFormat.Reader, AbstractString}, cachefile::AbstractString; kwargs...)
    if vcf isa AbstractString && isempty(basename(cachefile))
        cachefile = joinpath(cachefile, first(splitext(basename(vcf))) * ".jlso")
    end
    if isfile(cachefile)
        cdata = JLSO.load(cachefile, :vcfentries, :log)
        if !isempty(get(cdata, :log, ""))
            @warn "Found warnings in the log file, call `vcflog($(sprint(show, cachefile)))` for more information."
        end
        cdata[:vcfentries]
    else
        log = IOBuffer()
        logger = CompositeLogger(Logging.global_logger(), SimpleLogger(log, Warn))
        with_logger(logger) do
            vcfentries = IMPPROVE.vcfdata(vcf; kwargs...)
            JLSO.save(cachefile, :vcfentries => vcfentries,
                      :log => String(take!(log)))
            vcfentries
        end
    end
end

function vcflog(cachefile::AbstractString)
    JLSO.load(cachefile, :log)[:log]
end

end
