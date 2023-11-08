#!/usr/bin/env -S julia --startup-file=no
using Pkg
Pkg.activate(@__DIR__)

using CondaPkg

const SUPPORTED_HG_VERSIONS = ("19", "38")

function showhelp()
    println("""
    Liftover usage:
        $(basename(@__FILE__)) FROM TO INFILE OUTFILE REJECTFILE
    """)
end

if "--help" in ARGS
    showhelp()
    exit(0)
end

if length(ARGS) != 5
    println("! Wrong number of arguments, $(length(ARGS)) != 5")
    showhelp()
    exit(2)
end

hgfrom, hgto, infile, outfile, rejfile = ARGS

if hgfrom ∉ SUPPORTED_HG_VERSIONS
    println("! Unsupported Hg version: $hgfrom")
    showhelp()
    exit(3)
elseif hgto ∉ SUPPORTED_HG_VERSIONS
    println("! Unsupported Hg version: $hgfrom")
    showhelp()
    exit(3)
end

refseq(hgversion::String) =
    joinpath(@__DIR__, "data",
             Dict("19" => "GRCh37.p13.genome.fa.gz",
                  "38" => "GRCh38.primary_assembly.genome.fa.gz",
                  )[hgversion])

if endswith(infile, r"\.gvcf(?:\.gz)?")
    vcf_in_file = tempname(cleanup=false) * ".vcf"
    @info "Filtering out <NON_REF> terms in the GVCF"
    CondaPkg.withenv() do
        open(vcf_in_file, "w") do io
            for line in eachline(
                pipeline(`bcftools convert $infile
                          --fasta-ref $(abspath("./data/Hg38p14.fa")) --gvcf2vcf`,
                         `bcftools view -e 'ALT[0] == "<NON_REF>"'`))
                print(io, replace(line, ",<NON_REF>" => ""), '\n')
            end
        end
    end
    infile = vcf_in_file
end

CondaPkg.withenv() do
    chainfile = joinpath(@__DIR__, "data", "hg$(hgfrom)ToHg$(hgto).over.chain.gz")
    proc = `gatk LiftoverVcf
            --INPUT $infile --OUTPUT $outfile --REJECT $rejfile
            --CHAIN $chainfile --REFERENCE_SEQUENCE $(refseq(hgto))
            --LOG_FAILED_INTERVALS true --LIFTOVER_MIN_MATCH 1.0
            --WARN_ON_MISSING_CONTIG true` |>
            ignorestatus |> run
    exit(proc.exitcode)
end
