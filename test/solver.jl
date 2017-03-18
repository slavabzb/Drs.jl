#!/usr/bin/env julia

using MathProgBase
using Logging
using Drs

A = Float64[3 2 1; 2 5 3]
b = Float64[10, 15]
c = Float64[-2, -3, -4]

s = linprog(c, A, '<', b, -Inf, Inf, DrsMathProgSolver(LogLevel=DEBUG))

# The optimal solution value is Z = -20
# X1 = 0
# X2 = 0
# X3 = 5
