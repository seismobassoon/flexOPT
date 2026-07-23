#!/usr/bin/env julia

import Pkg

const FLEXOPT_ROOT = normpath(joinpath(@__DIR__, ".."))
const THREADS = get(ENV, "FLEXOPT_THREADS", "8")
const JULIA_MINOR = "$(VERSION.major).$(VERSION.minor)"
const SPECNAME = "flexopt-$(THREADS)-threads-$(JULIA_MINOR)"
const DISPLAYNAME = "flexOPT $(THREADS) threads $(JULIA_MINOR)"

VERSION >= v"1.12" && VERSION < v"1.13" ||
    error("flexOPT currently requires Julia 1.12.x; running $(VERSION)")

Pkg.activate(FLEXOPT_ROOT)
Pkg.instantiate()
Pkg.precompile()

import IJulia

kernel_path = IJulia.installkernel(
    "flexOPT $(THREADS) threads",
    "--startup-file=no",
    "--threads=$(THREADS)",
    "--project=$(FLEXOPT_ROOT)";
    specname=SPECNAME,
    displayname=DISPLAYNAME,
    env=Dict("JULIA_NUM_PRECOMPILE_TASKS" => "1"),
)

println("Installed Jupyter kernel:")
println("  $(DISPLAYNAME)")
println("  $(kernel_path)")
println("Project:")
println("  $(joinpath(FLEXOPT_ROOT, "Project.toml"))")
