const loadcache = Dict{Tuple, Any}()

function loadfromcache(name::AbstractString; kwargs...)
    key = (name, NamedTuple(kwargs))
    if !haskey(loadcache, key)
        loadcache[key] = IMPPROVE.getmach(name; kwargs...)
    end
    loadcache[key]
end

function loadfromcache(name::AbstractString, attrs::Symbol...; kwargs...)
    keys = tuple.(name, attrs, Ref(NamedTuple(kwargs)))
    if all(key -> haskey(loadcache, key), keys)
        getindex.(Ref(loadcache), keys)
    elseif haskey(loadcache, (name, NamedTuple(kwargs)))
        getindex.(Ref(loadcache[(name, NamedTuple(kwargs))]), attrs)
    else
        newkeys = filter(key -> !haskey(loadcache, key), keys)
        for (attr, val) in IMPPROVE.getmach(name, getindex.(newkeys, 2)...; kwargs...)
            loadcache[(name, attr, NamedTuple(kwargs))] = val
        end
        getindex.(Ref(loadcache), keys)
    end |> NamedTuple{attrs}
end

function plain!(axis::Axis)
    hidespines!(axis)
    hidedecorations!(axis)
end

function tookmsg(seconds::Number, prevline::Bool=true)
    if prevline
        print(stderr, "\e[F", "\e[", displaysize(stderr)[2] - 14 - (seconds >= 10), "C ")
    end
    println(stderr, "(took ", round(seconds, digits=2), "s)")
end

macro timed(msg, ex)
    quote
        try
        @info $(esc(msg))
        local starttime = time()
        local val = $(esc(ex))
        tookmsg(time() - starttime)
        val
        catch e
            @info "Error during: " * string($(esc(msg)))
            rethrow(e)
        end
    end
end

function Base.notify(@nospecialize(observable::AbstractObservable))
    val = observable[]
    for (_, f) in Observables.listeners(observable)::Vector{Pair{Int, Any}}
        result = try
            Base.invokelatest(f, val)
        catch err
            @info "Failure while trying to notify observable with " f observable err
            rethrow(err)
        end
        if result isa Consume && result.x
            # stop calling callbacks if event got consumed
            return true
        end
    end
    return false
end
