#!/usr/bin/env julia

using Drs
using Logging

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(DrsMathProgSolver(; LogLevel=DEBUG))

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(DrsMathProgSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(DrsMathProgSolver())
