#!/usr/bin/env -S julia --startup-file=no
using Pkg
Pkg.activate(@__DIR__)

using CondaPkg

const HEADER = "##INFO=<ID=NONCE,Number=1,Type=Integer,Description=\"Variant nonce (number only used once)\">"

function showhelp()
    println("""
    Annotate a VCF file with a nonce, to track variants across conersion.

    Usage:
        $(basename(@__FILE__)) FILE.vcf
    """)
end

if length(ARGS) < 1
    println("! Insufficent arguments, $(length(ARGS)) < 1")
end

file = first(ARGS)
isfile(file) || error("$file is not a file")

CondaPkg.withenv() do
    annotfile = tempname(cleanup=false) * ".tsv"
    open(annotfile, "w") do io
        for (i, line) in enumerate(eachline(`bcftools query -f "%CHROM\t%POS\n" $file`))
            write(io, line, '\t', string(i), '\n')
        end
    end
    run(`bgzip $annotfile`)
    annotfile = annotfile * ".gz"
    run(`tabix -s1 -b2 -e2 $annotfile`)
    run(`bcftools annotate
        -a $annotfile -H $HEADER
        -c CHROM,POS,INFO/NONCE
        $file
        -o $(get(ARGS, 2, "-"))`)
end
