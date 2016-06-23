# getting instance
LPSolver()

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

