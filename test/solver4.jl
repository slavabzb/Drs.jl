#!/usr/bin/env julia

# publication11

using Base.Test
using MathProgBase
using Logging
using Drs

A = Float64[1 1; 1 2; 3 1]
b = Float64[3, 5, 6]
c = Float64[-1, -2]

s = linprog(c, A, '<', b, -Inf, Inf, DrsMathProgSolver(LogLevel=DEBUG))

@test s.status == :Optimal
@test s.objval == -5
@test s.sol == [0, 5/2]
