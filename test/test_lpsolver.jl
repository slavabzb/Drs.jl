using MathProgBase
using LPLib

c = Float64[-3, -5]
b = Float64[180, 150, 300]
A = Float64[1 0; 0 2; 3 2]
sense = ['<', '<', '<']
lower_bounds = -Inf
upper_bounds = Inf

sol = linprog(c, A, sense, b, lower_bounds, upper_bounds, LPSolver())

@test sol.status == :Optimal
@test sol.objval == -525
@test sol.sol == [0 75 130 0 0]

println(sol)

exit()

A = Float64[0.4 0.2; 0.4 0.3; 0.2 0.5; 1.0 1.0]
b = Float64[10, 12, 10, 28]
c = Float64[-400, -300]

sense = ['<', '<', '<', '<']
lower_bounds = -Inf
upper_bounds = Inf

sol = linprog(c, A, sense, b, lower_bounds, upper_bounds, LPSolver(sense=:Max))

println(sol)

# testing equality
@test 1 == 1

# testing approximate equality
@test_approx_eq(1, 0.99999999999999)
@test_approx_eq_eps(1, 0.9, 0.1)

#testing throwing
function foo()
	error("error")
end

@test_throws ErrorException foo()

