module DrsMathProgSolverInterface

include("Simplex.jl")
using .Simplex

using Logging
@Logging.configure(level=DEBUG)

importall MathProgBase.SolverInterface

export DrsMathProgModel,
	DrsMathProgSolver,
	loadproblem!,
	optimize!,
	status,
	getreducedcosts,
	getconstrduals,
	getobjval,
	getsolution

immutable DrsMathProgSolver <: AbstractMathProgSolver
	options
end
DrsMathProgSolver(; kwargs...) = DrsMathProgSolver(kwargs)

type DrsMathProgModel <: AbstractLinearQuadraticModel
	A						# constraint coefficients
	b						# RHS
	c						# objective coefficients
	basis				# basis variables
	nonbasis		# nonbasis variables
end
LinearQuadraticModel(s::DrsMathProgSolver) = DrsMathProgModel(; s.options...)

function setparameters!(m::Union{DrsMathProgSolver, DrsMathProgModel}; kwargs...)
	for (option, value) in kwargs
		if option == :TimeLimit
			println("WARNING: TODO: $option")
		elseif option == :Silent && value == true
			Logging.configure(level=OFF)
		elseif option == :LogLevel
			Logging.configure(level=value)
		else
			println("WARNING: $option is unsupported")
		end
	end
end

function DrsMathProgModel(; kwargs...)
	m = DrsMathProgModel(0, 0, 0, 0, 0)
	setparameters!(m; kwargs...)
	m
end

function loadproblem!(m::DrsMathProgModel, A, l, u, c, lb, ub, sense)
	@debug("loadproblem!: A $A, l $l, u $u, c $c, lb $lb, ub $ub, sense $sense")
	m.A = A
	m.b = zeros(size(ub))
	m.c = c

	DrsTransformToStandardForm!(m, lb, ub, sense)

	r, c = size(m.A)
	m.basis = zeros(Int, r)

	DrsFindPotentialBasis!(m)
end

function DrsFindPotentialBasis!(m::DrsMathProgModel)
	r, c = size(m.A)
	for ic in 1:c
		column = m.A[:,ic]
		if countnz(column) == 1
			ir = findfirst(x -> x == 1, column)
			if ir != 0 && m.basis[ir] == 0
				# add the column if current row has not been selected
				m.basis[ir] = ic
			end
		end
	end
	m.nonbasis = setdiff(1:c, m.basis)
end

function DrsTransformToStandardForm!(m::DrsMathProgModel, lb, ub, sense)
	@assert length(lb) == length(ub) "the lengths of lower bounds ($(length(lb))) and upper bounds ($(length(ub))) are different"

	r, c = size(m.A)

	# check if b is negative

	for i in 1:length(lb)
		if lb[i] == typemin(typeof(lb[i]))
			# <
			if ub[i] < 0
				m.A[i,:] = -m.A[i,:]
				lb[i], ub[i] = -ub[i], -lb[i]
			end
		elseif ub[i] == typemax(typeof(ub[i]))
			# >
			if lb[i] < 0
				m.A[i,:] = -m.A[i,:]
				lb[i], ub[i] = -ub[i], -lb[i]
			end
		else
			# =
			if lb[i] < 0
				m.A[i,:] = -m.A[i,:]
				lb[i], ub[i] = -ub[i], -lb[i]
			end
		end
	end

	# add variables
	surplus = []

	for i in 1:length(lb)
		if lb[i] == typemin(typeof(lb[i]))
			# <, add slack
			m.A = [m.A zeros(r, 1)]
			m.A[i, end] = 1
			m.b[i] = ub[i]
			m.c = [m.c; 0]
		elseif ub[i] == typemax(typeof(ub[i]))
			# >, add surplus
			m.A = [m.A zeros(r, 1)]
			m.A[i, end] = -1
			m.b[i] = lb[i]
			m.c = [m.c; 0]
			push!(surplus, i)
		end
	end

	# =, add artificial
	for i in 1:length(lb)
		if lb[i] != typemin(typeof(lb[i])) && ub[i] != typemax(typeof(ub[i]))
			m.A = [m.A zeros(r, 1)]
			m.A[i, end] = 1
			m.b[i] = lb[i]
			m.c = [m.c; 0]
		end
	end

	# add artificial for surplus
	for i in surplus
		m.A = [m.A zeros(r, 1)]
		m.A[i, end] = 1
		m.c = [m.c; 0]
	end

	if sense == :Max
		m.c = -m.c
	end
end

function optimize!(m::DrsMathProgModel)
	@debug("A $(m.A)")
	@debug("b $(m.b)")
	@debug("c $(m.c)")
	@debug("basis $(m.basis)")
	@debug("nonbasis $(m.nonbasis)")

	B = m.A[:,m.basis]
	@debug("B $B")

	N = m.A[:,m.nonbasis]
	@debug("N $N")

	invB = SharedArray(typeof(B[1]), size(B),
		init = S -> S[linearindices(B)] = inv(B)[linearindices(B)])

	@debug("invB $invB")

	basic_vars = invB * m.b
	@debug("basic_vars $basic_vars")

	# find good candidates to leave the basis
	P = CHUZR(invB, basic_vars)
	@debug("P $P")

	pivotal_row = zeros(length(P))

	@sync begin
		for p in P
			pi = @spawn BTRAN(invB, p)
			pivotal_row = @spawn PRICE(N, fetch(pi))
			t = fetch(pivotal_row)
			@debug("t $t")
			# pivotal_row[p]
		end
	end

	@debug("pivotal_row $pivotal_row")
end

function status(m::DrsMathProgModel)
	@debug("status")
end

function getreducedcosts(m::DrsMathProgModel)
	@debug("getreducedcosts")
end

function getconstrduals(m::DrsMathProgModel)
	@debug("getconstrduals")
end

function getobjval(m::DrsMathProgModel)
	@debug("getobjval")
end

function getsolution(m::DrsMathProgModel)
	@debug("getsolution")
end

end
