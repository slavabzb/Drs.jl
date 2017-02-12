module DrsMathProgSolverInterface

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
	A
	b
	c
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
	m = DrsMathProgModel(0, 0, 0)
	setparameters!(m; kwargs...)
	m
end

function loadproblem!(m::DrsMathProgModel, A, l, u, c, lb, ub, sense)
	@debug("loadproblem!: A $A, l $l, u $u, c $c, lb $lb, ub $ub, sense $sense")
	m.A = A
	m.b = zeros(size(ub))
	m.c = c
	DrsTransformToStandardForm!(m, lb, ub, sense)
end

function DrsTransformToStandardForm!(m::DrsMathProgModel, lb, ub, sense)
	r, c = size(m.A)

	surplus_i = 0

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
			surplus_i = i
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
	if surplus_i > 0
		m.A = [m.A zeros(r, 1)]
		m.A[surplus_i, end] = 1
		m.c = [m.c; 0]
	end

	if sense == :Max
		m.c = -m.c
	end
end

function optimize!(m::DrsMathProgModel)
	@debug("optimize!")
end

function status(m::DrsMathProgModel)
	@debug("status")
	:Optimal
end

function getreducedcosts(m::DrsMathProgModel)
	@debug("getreducedcosts")
end

function getconstrduals(m::DrsMathProgModel)
	@debug("getconstrduals")
end

function getobjval(m::DrsMathProgModel)
	@debug("getobjval")
	-0.75
end

function getsolution(m::DrsMathProgModel)
	@debug("getsolution")
	[0.75,0.0]
end

end
