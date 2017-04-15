#!/usr/bin/env julia

# publication11

using Base.Test
using MathProgBase
using Logging
using Drs

A = Float64[3 4; 6 1]
b = Float64[6, 3]
c = Float64[-2, -1]

s = linprog(c, A, '<', b, -Inf, Inf, DrsMathProgSolver(LogLevel=DEBUG))
println(s)

@test s.status == :Optimal
@test s.objval == 13/7
@test s.sol == [2/7, 9/7]
