import Tables
using Base.Threads

struct MatrixTable{T}
    data::AbstractMatrix{T}
    names::Tuple{Vararg{Symbol}}
    lookup::Dict{Symbol, Int}
    lock::Threads.SpinLock
end

MatrixTable(data::Matrix{T}, names::Tuple{Vararg{Symbol}}) where {T} =
    MatrixTable{T}(data, names, Dict{Symbol, Int}(), Threads.SpinLock())
MatrixTable(data::Matrix, names::Tuple{Vararg{String}}) =
    MatrixTable(data, Symbol.(names))
MatrixTable(data::Matrix, names::AbstractVector{<:Union{Symbol, String}}) =
    MatrixTable(data, Tuple(names))
MatrixTable(data::Matrix) =
    MatrixTable(data, "x" .* string.(axes(data, 2)))

Base.names(mt::MatrixTable) = getfield(mt, :names)
Base.collect(mt::MatrixTable) = getfield(mt, :data)
function lookup(mt::MatrixTable, col::Symbol)
    if !haskey(getfield(mt, :lookup), col)
        lock(getfield(mt, :lock)) do
            getfield(mt, :lookup)[col] =
                @something(findfirst(==(col), names(mt)),
                           error("MatrixTable has no column '$col'"))
        end
    end
    getfield(mt, :lookup)[col]
end

Base.size(mt::MatrixTable) = size(collect(mt))
Base.size(mt::MatrixTable, I::Integer) = size(collect(mt), I)
Base.axes(mt::MatrixTable) = axes(collect(mt))
Base.axes(mt::MatrixTable, I::Integer) = axes(collect(mt), I)

for indexfn in (:getindex, :view)
    eval(quote
             Base.$(indexfn)(mt::MatrixTable, idx) =
                 $(indexfn)(collect(mt), idx)
             Base.$(indexfn)(mt::MatrixTable, row::Union{<:Integer, <:AbstractVector{<:Integer}, Colon},
                           col::Integer) =
                 $(indexfn)(collect(mt), row, col)
             Base.$(indexfn)(mt::MatrixTable, row::Union{<:Integer, <:AbstractVector{<:Integer}, Colon},
                           col::Colon) =
                 MatrixTable($(indexfn)(collect(mt), row, col),
                             getfield(mt, :names),
                             getfield(mt, :lookup),
                             getfield(mt, :lock))
             Base.$(indexfn)(mt::MatrixTable, row::Union{<:Integer, <:AbstractVector{<:Integer}, Colon},
                           col::AbstractVector{<:Integer}) =
                 MatrixTable($(indexfn)(collect(mt), row, col),
                             getfield(mt, :names)[col],
                             Dict{Symbol, Int}(), Threads.SpinLock())
             Base.$(indexfn)(mt::MatrixTable, row::Any, col::Symbol) =
                 $(indexfn)(mt, row, lookup(mt, col))
             Base.$(indexfn)(mt::MatrixTable, row::Any, cols::Vector{Symbol}) =
                 $(indexfn)(mt, row, lookup.(Ref(mt), cols))
             Base.$(indexfn)(mt::MatrixTable, row::Any, cols::Vector{<:AbstractString}) =
                 getindex(mt, row, lookup.(Ref(mt), Symbol.(cols)))
         end)
end

Base.getproperty(mt::MatrixTable, col::Symbol) =
    collect(mt)[:, lookup(mt, col)]

Tables.istable(::Type{<:MatrixTable}) = true

Tables.schema(mt::MatrixTable{T}) where {T} =
    Tables.Schema{nothing, nothing}(
        collect(names(mt)), fill(T, size(mt, 2)))

Tables.columnaccess(::Type{<:MatrixTable}) = true
Tables.columns(mt::MatrixTable) = mt

Tables.getcolumn(mt::MatrixTable, ::Type, col::Int, _::Symbol) = mt[:, col]
Tables.getcolumn(mt::MatrixTable, name::Symbol) = mt[:, name]
Tables.getcolumn(mt::MatrixTable, i::Int) = mt[:, i]
Tables.columnnames(mt::MatrixTable) = names(mt)

Tables.rowaccess(::Type{<:MatrixTable}) = true
Tables.rows(mt::MatrixTable) = mt

struct MatrixRow{T} <: Tables.AbstractRow
    row::Int
    source::MatrixTable{T}
end

Base.getindex(mr::MatrixRow, col::Union{<:Integer, <:AbstractVector{<:Integer}, Colon, Symbol}) =
    getfield(mr, :source)[getfield(mr, :row), col]
Base.getproperty(mr::MatrixRow, col::Symbol) =
    getindex(mr, col)

Base.eltype(::MatrixTable{T}) where {T} = MatrixRow{T}
Base.length(mt::MatrixTable) = size(mt, 1)

Tables.getcolumn(mr::MatrixRow, ::Type, col::Int, ::Symbol) =
    mr[col]
Tables.getcolumn(mr::MatrixRow, col::Int) = mr[col]
Tables.getcolumn(mr::MatrixRow, col::Symbol) = mr[col]
Tables.columnnames(mr::MatrixRow) = names(getfield(mr, :source))

Base.Matrix(mt::MatrixTable) = collect(mt)

MatrixTable(df::DataFrame) =
    MatrixTable(Matrix(df), names(df))

MMI.matrix(mt::MatrixTable) = collect(mt)
