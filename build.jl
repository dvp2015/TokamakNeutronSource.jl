#!/bin/bash
#=
exec julia --color=yes --startup-file=no -e 'include(popfirst!(ARGS))' "${BASH_SOURCE[0]}" "$@"
=#

using Pkg
Pkg.activate(".")
Pkg.build(; verbose = true)
Pkg.test(coverage=true)
