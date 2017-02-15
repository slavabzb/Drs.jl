#!/usr/bin/env julia

using Base.Test
using Drs.DrsMathProgSolverInterface

A = [2 3 1 0 0 0; -3 2 0 1 0 0; 0 2 0 0 1 0; 2 1 0 0 0 1]
b = nothing
c = nothing

row, col = size(A)

basis = zeros(Int, row)
nonbasis = nothing

m = DrsMathProgModel(A, b, c, basis, nonbasis)
Drs.DrsMathProgSolverInterface.DrsFindPotentialBasis!(m)

@test m.basis == [3, 4, 5, 6]
@test m.nonbasis == [1, 2]
