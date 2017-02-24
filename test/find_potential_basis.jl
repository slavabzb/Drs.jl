#!/usr/bin/env julia

using Base.Test
using Drs.DrsMathProgSolverInterface

A = [2 3 1 0 0 0; -3 2 0 1 0 0; 0 2 0 0 1 0; 2 1 0 0 0 1]

row, col = size(A)

m = DrsMathProgModel(A, nothing, nothing, zeros(Int, row), nothing, nothing, nothing, nothing)
Drs.DrsMathProgSolverInterface.DrsFindPotentialBasis!(m)

@test m.basis == [3, 4, 5, 6]
@test m.nonbasis == [1, 2]