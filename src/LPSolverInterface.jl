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
	A			# Prob data: coefficient matrix for constraints
	b			# Prob data: rhs
	c			# Prob data: coefficient for objective function
	
	inq			# m dimentional vector; 1 for >; 0 for =; -1 for <
	
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
	terminate	# Flag for termination: true - terminate; false - continue
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
	println("A ", A)
	println("collb ", collb)
	println("colub ", colub)
	println("obj ", obj)
	println("rowlb ", rowlb)
	println("rowub ", rowub)
	println("sense ", sense)
	println("model ", m)
	
	m.sense = sense
	
	prepare!(m, rowlb, rowub)
	
	m.A = A
	println("m.A ", m.A)
	
	checkNegative!(m)
	
	m.m, m.n = size(m.A)
	println("m.m ", m.m)
	println("m.n ", m.n)
	
	m.c = obj
	println("m.c ", m.c)
	
	m.M = ceil(max(norm(m.A), norm(m.b))/100)*100
	println("m.M ", m.M)
	
	m.basis = zeros(1, m.m)
	println("m.basis ", m.basis)
	
	m.terminate = false
	m.counter = 0
	
	transformToStandardForm!(m)
	initialize!(m)
	
	println("m.b ", m.b)
end

function prepare!(m::LPMathProgModel, rowlb, rowub)
	# Merge lower_bounds and upper_bounds to single b
	# Set inq array
	
	n = length(rowlb)
	
	m.inq = zeros(n)
	m.b = rowub
	
	for i in 1:n
		if rowlb[i] == -Inf
			m.inq[i] = -1
		elseif rowub[i] == Inf
			m.inq[i] = 1
			m.b[i] = rowlb[i]
		end
	end

	println("m.inq ", m.inq)
end

function checkNegative!(m::LPMathProgModel)
	# Check if b is negative
	
	b_Neg = [x ? 1 : 0 for x in (m.b .< 0)]
	println("b_Neg ", b_Neg)
	
	for (i, x) in enumerate(b_Neg)
		if x == 1
			m.inq[i] = -m.inq[i]
			m.b[i] = -m.b[i]
			m.A[i,:] = -m.A[i,:]
		end
	end
end

function transformToStandardForm!(m::LPMathProgModel)
	# Add slack variables
	
	n_slack = sum(m.inq .< 0) + sum(m.inq .> 0)
	println("n_slack ", n_slack)
	
	m.A = [m.A zeros(m.m, n_slack)]
	println("m.A ", m.A)
	
	m.c = [m.c; zeros(n_slack, 1)]
	println("m.c ", m.c)
	
	idx_slack = m.n + 1
	for i in 1:m.m
		if m.inq[i] == 1
			m.A[i, idx_slack] = -1
		elseif m.inq[i] == -1
			m.A[i, idx_slack] = 1
		end
		idx_slack += 1
	end
	
	updateDimension!(m)
end

function updateDimension!(m::LPMathProgModel)
	m.n = size(m.A)[2]
	println("updateDimension m.n ", m.n)
end

function initialize!(m::LPMathProgModel)
	# Find potential basis
	for i in 1:m.n
		column = m.A[:,i]
		if countnz(column) == 1
			row_number = find(x -> x == 1, m.A[:,i])
			if !isempty(row_number)
				# Add the col if the current row has not been selected
				if m.basis[row_number] == [0]
					m.basis[row_number] = i;
				end
			end
		end
	end
	println("m.basis ", m.basis)
	
	# Add big M if necessary
	n_artfVar = m.m - sum(m.basis .> 0)
	println("n_artfVar ", n_artfVar)
	
	if n_artfVar > 0
		# Record the index for artificial variables
		m.A = [m.A zeros(m.m, n_artfVar)]
		
		add_to_rows = find(x -> x == 0, m.basis)
		
		# Formulate matrix A with newly added artificial vars
		for i in 1:length(add_to_rows)
			m.A[add_to_rows[i], m.n + i] = 1
			m.basis[add_to_rows[i]] = m.n + i
		end
		
		m.c = [m.c; m.M * ones(n_artfVar, 1)]
	end
	
	m.nonbasis = setdiff(1:m.m, m.basis)
	updateDimension!(m)
	
	m.x = zeros(m.n, 1)
	println("m.x ", m.x)
	
	basis_ids = trunc(Int, m.basis[:])
	println("basis_ids ", basis_ids)
	
	println("lhs ", m.A[:,basis_ids])
	println("rhs ", m.b)
	
	m.x[basis_ids] = m.A[:,basis_ids] \ m.b
	println("m.x ", m.x)
	
	m.d = m.c - m.A' * m.c[basis_ids]
	println("m.d ", m.d)
	
	m.z = 0
end

function incrementCounter!(m::LPMathProgModel)
	m.counter += 1
	if m.counter == m.maxIter
		m.terminate = true
		m.status = :UserLimit
	end
end

function computeDualVars!(m::LPMathProgModel)
	# BTRAN
	basis_ids = trunc(Int, m.basis[:])
	m.y = inv(m.A[:,basis_ids])' * m.c[basis_ids]
	println("m.y ", m.y)
end

function priceNonBasicVars!(m::LPMathProgModel)
	# PRICE
	basis_ids = trunc(Int, m.basis[:])
	nonbasis_ids = trunc(Int, m.nonbasis[:])
	m.d[nonbasis_ids] = m.c[nonbasis_ids] - m.A[:,nonbasis_ids]' * m.y
	m.d[basis_ids] = 0
	println("m.d ", m.d)
end

function chooseNonBasisVarToEnterBasis!(m::LPMathProgModel)
	# CHUZC
	nonbasis_ids = trunc(Int, m.nonbasis[:])
	m.dNq = minimum(m.d[nonbasis_ids])
	if m.dNq >= 0
		m.terminate = true
		m.status = :Optimal
		if m.sense == :Max
			m.fval = -m.z
		else
			m.fval = m.z
		end
	else
		m.q = findfirst(x -> x == m.dNq, m.d[nonbasis_ids])
		m.q = Int(m.nonbasis[m.q])
		println("m.q ", m.q)
	end
end

function findUpdatedCol!(m::LPMathProgModel)
	# FTRAN
	basis_ids = trunc(Int, m.basis[:])
	m.updCol = m.A[:,basis_ids] \ m.A[:,m.q]
	println("m.updCol ", m.updCol)
end

function findBasisVarToLeaveBasis!(m::LPMathProgModel)
	# CHUZR
	basis_ids = trunc(Int, m.basis[:])
	m.sigma, indx = findmin(m.x[basis_ids] ./ m.updCol)
	println("m.sigma ", m.sigma)
	println("indx ", indx)
	
	if isinf(m.sigma)
		m.terminate = true
		m.status = :Unbounded
	end
	
	m.t = m.basis[indx]
	println("m.t ", m.t)
end

function updateStep!(m::LPMathProgModel)
	e = zeros(m.n, 1)
	e[m.q] = 1
	
	basis_ids = trunc(Int, m.basis[:])
	nonbasis_ids = trunc(Int, m.nonbasis[:])
	
	m.x[basis_ids] = m.x[basis_ids] - m.sigma * m.updCol
	m.x[nonbasis_ids] = m.x[nonbasis_ids] + m.sigma * e[nonbasis_ids]
	
	m.z = m.z + m.dNq * m.sigma
	println("m.z ", m.z)
end

function updateBasis!(m::LPMathProgModel)
	m.basis[find(x -> x == m.t, m.basis)] = m.q
	println("m.basis ", m.basis)
	m.nonbasis[find(x -> x == m.q, m.nonbasis)] = m.t
	println("m.nonbasis ", m.nonbasis)
end

function optimize!(m::LPMathProgModel)
	println("optimize ")
	
	while !m.terminate
		incrementCounter!(m)
		computeDualVars!(m)
		priceNonBasicVars!(m)
		chooseNonBasisVarToEnterBasis!(m)
		findUpdatedCol!(m)
		findBasisVarToLeaveBasis!(m)
		updateStep!(m)
		updateBasis!(m)
	end
end

function status(m::LPMathProgModel)
	println("status ", m.status)
	println("m.x ", m.x)
	println("m.fval ", m.fval)
end

end

