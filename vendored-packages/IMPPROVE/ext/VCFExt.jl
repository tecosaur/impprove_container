module VCFExt

using IMPPROVE

using VariantCallFormat
using InlineStrings
using DataFrames
using TranscodingStreams
using CodecZlib

import IMPPROVE.vcfdata

"""
    vcfdata(vcf::VariantCallFormat.Reader, annotator::Function=identity)
Extract all records from `vcf`, annotating with extra information by calling
`annotator` on the DataFrame produced.
"""
function vcfdata(vcf::VariantCallFormat.Reader; annotator::Function=identity, filter::Function=Returns(true))
    vcfd = DataFrame(chromosome=String3[], location=Int[], change=Any[])
    ChromStr = String3
    for record in vcf
        if all(!isempty, (record.ref, record.alt)) && filter(record)
            ref = VariantCallFormat.ref(record)
            alt = VariantCallFormat.alt(record)
            chrom = replace(VariantCallFormat.chrom(record), "chr" => "")
            if ChromStr == String3 && ncodeunits(chrom) > 3
                @warn "Unexpectadly long chromosome name seen: $(VariantCallFormat.chrom(record))"
                ChromStr = String
                vcfd.chromosome = String.(vcfd.chromosome)
            end
            if all(==(1), length.((ref, alt)))
                push!(vcfd,
                    (chromosome = ChromStr(chrom),
                     location = VariantCallFormat.pos(record),
                     change = Pair(ref[1], alt[1][1])))
            end
        end
    end
    annotator(vcfd)
end

"""
    vcfdata(f::String, args...; kargs...)
Open `f` with `VariantCallFormat.Reader`, and call `vcfdata` with it and all
provided arguments.
"""
vcfdata(f::String, args...; kargs...) =
    vcfdata(VariantCallFormat.Reader(
        if endswith(f, ".gz")
            TranscodingStream(GzipDecompressor(), open(f, "r"))
        else
            open(f, "r")
        end), args...; kargs...)

end
