module LPSolverInterface

importall MathProgBase.SolverInterface

export LPSolver
immutable LPSolver <: AbstractMathProgSolver
	options
end
LPSolver(;kwargs...) = LPSolver(kwargs)

type LPMathProgModel <: AbstractLinearQuadraticModel
	options
	
	# Problem data
	A			# Porb data: coefficient matrix for constraints
	b			# Porb data: rhs
	c			# Porb data: coefficient for objective function
	
	inq			# m dimentional vector; 1 for >=; 0 for ==; -1 for <=
	
	m			# Number of constraints
	n			# Number of variables
	fval		# Optimal objective value
	sense		# Min or max
	
	# Iteration info
	x			# Vector of primal variables
	y			# Vector of dual variables
	z			# Current objective value
	
	d			# Vector of reduced costs
	dNq			# The value of the most negative reduced cost
	q			# Nonbasic variable to enter the basis
	t			# Basic variable to leave the basis
	
	updCol		# Updated column
	
	sigma		# Step size
	basis		# Current basis
	nonbasis	# Current nonbasis
	M			# Big M
	counter		# Iteration counter
	terminate	# Flag for termination: 1 terminate; 0 continue
	status		# maxIter, unbounded, infeasible, optimal
	
	maxIter		# Max number of iterations
	
	function LPMathProgModel(;options...)
		model = new()
		model.options = options
		model.maxIter = 20
		model
	end
end
LinearQuadraticModel(s::LPSolver) = LPMathProgModel(;s.options...)

function loadproblem!(m::LPMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
	println("loadproblem")
	
	m.sense = sense
	println("m.sense ", m.sense)
	
	if m.sense == Symbol("Max")
		obj = -obj
		println("obj ", obj)
	end
	
	b_Neg = rowub .< 0
	println("b_Neg ", b_Neg)
	
	m.A = A
	println("m.A ", m.A)
	
	m.b = rowub
	println("m.b ", m.b)
	
	m.c = obj
	println("m.c ", m.c)
	
	m.M = ceil(max(norm(m.A), norm(m.b))/100)*100
	println("m.M ", m.M)
	
	m.m, m.n = size(m.A)
	println("m.m ", m.m)
	println("m.n ", m.n)
	
	m.basis = zeros(1, m.m)
	println("m.basis ", m.basis)
	
	m.terminate = 0
	m.counter = 0
end

function optimize!(m::LPMathProgModel)
println("optimize")
end

function status(m::LPMathProgModel)
println("status")
end

end

