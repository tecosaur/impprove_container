function (; path::String)
    @import HDF5: h5open, readmmap, ismmappable
    @import Statistics: mean
    cadd = h5open(path, "r")
    function getcadd(chromosome::AbstractString, location::Integer,
                        offset::Integer, score::AbstractString = "raw")
        index = searchsortedfirst(axes(cadd[chromosome]["positions"], 1), (),
                                    by = x -> if x isa Int
                                        cadd[chromosome]["positions"][x, 1]
                                    else location end)
        cadd[chromosome][score][index+offset,1]
    end
    base_order = ['A', 'C', 'G', 'T']
    multi_bases = Dict(
        'R' => ['G', 'A'],
        'Y' => ['C', 'T'],
        'M' => ['A', 'C'],
        'K' => ['G', 'T'],
        'S' => ['G', 'C'],
        'W' => ['A', 'T'],
        'B' => ['C', 'G', 'T'],
        'D' => ['A', 'G', 'T'],
        'H' => ['A', 'C', 'T'],
        'V' => ['A', 'C', 'G'],
        'N' => ['A', 'C', 'G', 'T'])
    change_offset_map = Dict(
        (b2 => b1) => i - 1
        for b2 in base_order
            for (i, b1) in enumerate(setdiff(base_order, b2)))
    function getcadd(chromosome::AbstractString, location::Integer,
                     change::Pair{Char,Char}, score::AbstractString = "raw")
        if first(change) == last(change)
            0.0
        elseif last(change) ∈ ('A', 'C', 'G', 'T')
            getcadd(chromosome, location, change_offset_map[change], score)
        elseif haskey(multi_bases, last(change))
            getcadd.(chromosome, location,
                     first(change) .=> multi_bases[last(change)], score) |> mean
        else
            throw(KeyError(change))
        end
    end
    function getcadd(chromosome::AbstractString, location::Integer,
                     score::AbstractString = "raw")
        getcadd.(chromosome, location, (0, 1, 2), score) |> mean
    end
    getcadd
end
