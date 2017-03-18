#!/usr/bin/env julia

using Logging
using Base.Test
using Drs.DrsMathProgSolverInterface

Logging.configure(level=DEBUG)

A = [2 3 1 0 0 0; -3 2 0 1 0 0; 0 2 0 0 1 0; 2 1 0 0 0 1]
b = [6, 3, 5, 4]
c = nothing

basis = [3, 4, 5, 6]
nonbasis = [1, 2]

m = DrsMathProgModel(A, b, c, basis, nonbasis, nothing, nothing, nothing)
Drs.DrsMathProgSolverInterface.DrsChuzr!(m)
