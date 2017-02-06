module DrsMathProgSolverInterface

using Logging
@Logging.configure(level=DEBUG)

include("DrsCore.jl")
using .DrsCore

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
	stub
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
	m = DrsMathProgModel(0)
	setparameters!(m; kwargs...)
	m
end

function loadproblem!(m::DrsMathProgModel, A, l, u, c, lb, ub, sense)
	@debug("loadproblem!: A $A, l $l, u $u, c $c, lb $lb, ub $ub, sense $sense")
end

function optimize!(m::DrsMathProgModel)
	@debug("optimize!")
	CHUZR()
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
