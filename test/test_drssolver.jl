#!/usr/bin/env julia

using Base.Test
using MathProgBase
using Logging
using Drs

logLevel=OFF

A = Float64[1 0; 0 2; 3 2]
b = Float64[180, 150, 300]
c = Float64[-3, -5]

sense = ['<', '<', '<']
lower_bounds = -Inf
upper_bounds = Inf

sol = linprog(c, A, sense, b, lower_bounds, upper_bounds, DrsSolver(logLevel=logLevel))

@test sol.status == :Optimal
@test sol.objval == -525
@test sol.sol == [0 75 130 0 0]



A = Float64[0.4 0.2; 0.4 0.3; 0.2 0.5; 1.0 1.0]
b = Float64[10, 12, 10, 28]
c = Float64[-400, -300]

sense = ['<', '<', '<', '<']
lower_bounds = -Inf
upper_bounds = Inf

sol = linprog(c, A, sense, b, lower_bounds, upper_bounds, DrsSolver(logLevel=logLevel))

@test sol.status == :Optimal
@test sol.objval == -10600
@test_approx_eq_eps(sol.sol, [22 0 0 1.4 2.6 0], 1e-6)



A = Float64[1 4; 2 3; 2 1]
b = Float64[4, 6, 4]
c = Float64[-2, -3]

sense = ['<', '<', '<']
lower_bounds = -Inf
upper_bounds = Inf

sol = linprog(c, A, sense, b, lower_bounds, upper_bounds, DrsSolver(logLevel=logLevel))

@test sol.status == :Optimal
@test_approx_eq(sol.objval, -5.142857142857142)
@test_approx_eq_eps(sol.sol, [0 0.571429 0 0.857143 0], 1e-6)

