#!/usr/bin/env -S julia --startup-file=no
using Pkg
Pkg.activate(@__DIR__)

using CondaPkg

CondaPkg.withenv() do
    proc = vcat("gatk", ARGS) |> Cmd |> ignorestatus |> run
    exit(proc.exitcode)
end
