module GeneVariantsMinimalData

using DataToolkit
using DataFrames, CSV

function __init__()
    DataToolkit.init(@__MODULE__, force=true)
    DataToolkit.@addpkgs CodecZstd CSV HDF5 Statistics
end

# Util functions

function find_surrounding_genes(location::UnitRange, chromosome, table::DataFrame,
                                locationkey::Symbol=:location, chromosomekey::Symbol=:chromosome)
    table_chromosone = filter(chromosomekey => c -> c == chromosome, table)
    pos = findall(l -> location ⊆ l, table_chromosone[!, locationkey])
    table_chromosone[pos, :]
end

find_surrounding_genes(location::Integer, chromosome, table::DataFrame,
                       locationkey::Symbol=:location, chromosomekey::Symbol=:chromosome) =
                           find_surrounding_genes(location:location, chromosome, table, locationkey, chromosomekey)

function find_surrounding_exons_genes(location::Union{Integer, UnitRange}, chromosome, exon_table::DataFrame,
                                      locationkey::Symbol=:location, chromosomekey::Symbol=:chromosome;
                                      protprefer::Bool=false, protfilter::Bool=false)
    exons = find_surrounding_genes(
        location, chromosome, exon_table, locationkey, chromosomekey)
    if protfilter || protprefer && "protein_coding" in getindex.(exons.attributes, "gene_type")
        filter!(:attributes => a -> a["gene_type"][1] == "protein_coding",
                exons)
    end
    map(e -> replace(e["gene_id"][1], r"\..*$" => ""), exons.attributes) |>
        unique
end

end
