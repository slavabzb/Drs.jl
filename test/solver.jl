#!/usr/bin/env julia

using MathProgBase
using Logging
using Drs

# A = Float64[1 0; 0 2; 3 2]
# b = Float64[180, 150, 300]
# c = Float64[-3, -5]

A = Float64[3 5 2; 4 4 4; 2 4 5]
b = Float64[60, 72, 100]
c = Float64[5, 10, 8]

s = linprog(c, A, ['=', '>', '<'], b, -Inf, Inf, DrsMathProgSolver(LogLevel=DEBUG))
