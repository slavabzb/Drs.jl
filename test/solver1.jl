#!/usr/bin/env julia

using Base.Test
using MathProgBase
using Logging
using Drs

A = Float64[3 2 1; 2 5 3]
b = Float64[10, 15]
c = Float64[-2, -3, -4]

s = linprog(c, A, '<', b, -Inf, Inf, DrsMathProgSolver())

@test s.status == :Optimal
@test s.objval == -20
@test s.sol == [0, 0, 5]
