#!/usr/bin/env julia

using Base.Test
using Drs.DrsMathProgSolverInterface

# Test case
A = [1 0; 0 2; 3 2]
b = [180, 150, 300]
c = [-3, -5]

lb = [-Inf, 150, 300]
ub = [180, 150, Inf]

mA = [1 0 1 0 0 0;
      0 2 0 0 1 0;
      3 2 0 -1 0 1]
mc = [-3, -5, 0, 0, 0, 0]

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Min)

@test m.A == mA
@test m.b == b
@test m.c == mc

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Max)

@test m.A == mA
@test m.b == b
@test m.c == -mc


# Test case
A = [1 0; 0 2; 3 2]
b = [180, 150, 300]
c = [-3, -5]

lb = [180, 150, 300]
ub = [Inf, Inf, Inf]

mA = [1 0 -1 0 0 1 0 0;
      0 2 0 -1 0 0 1 0;
      3 2 0 0 -1 0 0 1]
mc = [-3, -5, 0, 0, 0, 0, 0, 0]

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Min)

@test m.A == mA
@test m.b == b
@test m.c == mc

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Max)

@test m.A == mA
@test m.b == b
@test m.c == -mc


# Test case
A = [1 0; 0 2; 3 2]
b = [-180, -150, -300]
c = [-3, -5]

lb = [-Inf, -150, -300]
ub = [-180, -150, Inf]

mA = [-1 0 -1 0 0 1;
      0 -2 0 0 1 0;
      -3 -2 0 1 0 0]
mb = [180, 150, 300]
mc = [-3, -5, 0, 0, 0, 0]

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Min)

@test m.A == mA
@test m.b == mb
@test m.c == mc


# Test case
A = [1 0; 0 2; 3 2]
b = [0, 0, 0]
c = [-3, -5]

lb = [-Inf, 0, 0]
ub = [0, 0, Inf]

mA = [1 0 1 0 0 0;
      0 2 0 0 1 0;
      3 2 0 -1 0 1]
mb = [0, 0, 0]
mc = [-3, -5, 0, 0, 0, 0]

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Min)

@test m.A == mA
@test m.b == mb
@test m.c == mc


# Test case
A = [3 5 2; 4 4 4; 2 4 5]
b = [60, 72, 100]
c = [5, 10, 8]

lb = [60, 72, -Inf]
ub = [60, Inf, 100]

mA = [3 5 2 0 0 1 0; 4 4 4 -1 0 0 1; 2 4 5 0 1 0 0]
mb = [60, 72, 100]
mc = [5, 10, 8, 0, 0, 0, 0]

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Min)

@test m.A == mA
@test m.b == mb
@test m.c == mc

m = DrsMathProgModel(A, b, c)
Drs.DrsMathProgSolverInterface.DrsTransformToStandardForm!(m, lb, ub, :Max)

@test m.A == mA
@test m.b == mb
@test m.c == -mc
