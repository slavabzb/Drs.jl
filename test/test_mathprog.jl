#!/usr/bin/env julia

using Drs

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(DrsSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(DrsSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(DrsSolver())
