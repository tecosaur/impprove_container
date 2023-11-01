#!/usr/bin/env julia -S --project=/
using bcftools_jll

open("bcftools", "w") do io
    println(io, "#!/bin/sh")
    cmd = bcftools()
    print(io, join(cmd.env, ' '))
    print(io, " exec ", join(cmd.exec, ' '))
    println(io, " \"\$@\"")
end
