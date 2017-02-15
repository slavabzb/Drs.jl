#!/usr/bin/env julia

using MathProgBase
using Logging
using Drs

A = Float64[2 1; 2 3; 3 1]
b = Float64[18, 42, 24]
c = Float64[3, 2]

s = linprog(c, A, ['<', '<', '<'], b, -Inf, Inf, DrsMathProgSolver(LogLevel=DEBUG))
