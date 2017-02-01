module LPSolverInterface

using Logging
@Logging.configure(level=DEBUG)

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
		setparams!(model)
		model
	end
	
	function setparams!(m::LPMathProgModel)
		for (name,value) in m.options
			if name == :logLevel
				Logging.configure(level=value)
			end
		end
	end
end

LinearQuadraticModel(s::LPSolver) = LPMathProgModel(;s.options...)

function loadproblem!(m::LPMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
	@debug("loadproblem")
	@debug("A ", A)
	@debug("collb ", collb)
	@debug("colub ", colub)
	@debug("obj ", obj)
	@debug("rowlb ", rowlb)
	@debug("rowub ", rowub)
	@debug("sense ", sense)
	@debug("model ", m)
	
	m.sense = sense
	
	prepare!(m, rowlb, rowub)
	
	m.A = A
	@debug("m.A ", m.A)
	
	checkNegative!(m)
	
	m.m, m.n = size(m.A)
	@debug("m.m ", m.m)
	@debug("m.n ", m.n)
	
	m.c = obj
	@debug("m.c ", m.c)
	
	m.M = ceil(max(norm(m.A), norm(m.b))/100)*100
	@debug("m.M ", m.M)
	
	m.basis = zeros(Int, m.m)
	@debug("m.basis ", m.basis)
	
	m.terminate = false
	m.counter = 0
	
	transformToStandardForm!(m)
	initialize!(m)
	
	@debug("m.b ", m.b)
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

	@debug("m.inq ", m.inq)
end

function checkNegative!(m::LPMathProgModel)
	# Check if b is negative
	
	b_Neg = [x ? 1 : 0 for x in (m.b .< 0)]
	@debug("b_Neg ", b_Neg)
	
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
	@debug("n_slack ", n_slack)
	
	m.A = [m.A zeros(m.m, n_slack)]
	@debug("m.A ", m.A)
	
	m.c = [m.c; zeros(n_slack, 1)]
	@debug("m.c ", m.c)
	
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
	@debug("updateDimension m.n ", m.n)
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
	@debug("m.basis ", m.basis)
	
	# Add big M if necessary
	n_artfVar = m.m - sum(m.basis .> 0)
	@debug("n_artfVar ", n_artfVar)
	
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
	@debug("m.x ", m.x)
	
	@debug("lhs ", m.A[:,m.basis])
	@debug("rhs ", m.b)
	
	m.x[m.basis] = m.A[:,m.basis] \ m.b
	@debug("m.x ", m.x)
	
	m.d = m.c - m.A' * m.c[m.basis]
	@debug("m.d ", m.d)
	
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
	m.y = inv(m.A[:,m.basis])' * m.c[m.basis]
	@debug("m.y ", m.y)
end

function priceNonBasicVars!(m::LPMathProgModel)
	# PRICE
	m.d[m.nonbasis[:]] = m.c[m.nonbasis[:]] - m.A[:,m.nonbasis[:]]' * m.y
	m.d[m.basis] = 0
	@debug("m.d ", m.d)
end

function chooseNonBasisVarToEnterBasis!(m::LPMathProgModel)
	# CHUZC
	m.dNq = minimum(m.d[m.nonbasis[:]])
	if m.dNq >= 0
		m.terminate = true
		m.status = :Optimal
		if m.sense == :Max
			m.fval = -m.z
		else
			m.fval = m.z
		end
	else
		m.q = findfirst(x -> x == m.dNq, m.d[m.nonbasis[:]])
		m.q = m.nonbasis[m.q]
		@debug("m.q ", m.q)
	end
end

function findUpdatedCol!(m::LPMathProgModel)
	# FTRAN
	m.updCol = m.A[:,m.basis] \ m.A[:,m.q]
	@debug("m.updCol ", m.updCol)
end

function findBasisVarToLeaveBasis!(m::LPMathProgModel)
	# CHUZR
	m.sigma, indx = findmin(m.x[m.basis] ./ m.updCol)
	@debug("m.sigma ", m.sigma)
	@debug("indx ", indx)
	
	if isinf(m.sigma)
		m.terminate = true
		m.status = :Unbounded
	end
	
	m.t = m.basis[indx]
	@debug("m.t ", m.t)
end

function updateStep!(m::LPMathProgModel)
	e = zeros(m.n, 1)
	e[m.q] = 1
	
	m.x[m.basis] = m.x[m.basis] - m.sigma * m.updCol
	m.x[m.nonbasis[:]] = m.x[m.nonbasis[:]] + m.sigma * e[m.nonbasis[:]]
	
	m.z = m.z + m.dNq * m.sigma
	@debug("m.z ", m.z)
end

function updateBasis!(m::LPMathProgModel)
	m.basis[find(x -> x == m.t, m.basis)] = m.q
	@debug("m.basis ", m.basis)
	m.nonbasis[find(x -> x == m.q, m.nonbasis)] = m.t
	@debug("m.nonbasis ", m.nonbasis)
end

function optimize!(m::LPMathProgModel)
	@debug("optimize ")
	
	while !m.terminate
		incrementCounter!(m)
		computeDualVars!(m)	# BTRAN
		priceNonBasicVars!(m)	# PRICE
		chooseNonBasisVarToEnterBasis!(m)	# CHUZC
		findUpdatedCol!(m)	# FTRAN
		findBasisVarToLeaveBasis!(m)	# CHUZR
		updateStep!(m)
		updateBasis!(m)
	end
end

function status(m::LPMathProgModel)
	@debug("status ", m.status)
	@debug("m.x ", m.x)
	@debug("m.fval ", m.fval)
	m.status
end

function getreducedcosts(m::LPMathProgModel)
end

function getconstrduals(m::LPMathProgModel)
end

function getobjval(m::LPMathProgModel)
	m.fval
end

function getsolution(m::LPMathProgModel)
	m.x'
end

end

