module LPSolverInterface

importall MathProgBase.SolverInterface

export LPSolver

type LPMathProgModel <: AbstractLinearQuadraticModel
	dummy
end

immutable LPSolver <: AbstractMathProgSolver
    options
end
LPSolver(;kwargs...) = LPSolver(kwargs)

function LPMathProgModel(;options...)
	m::LPMathProgModel = LPMathProgModel(0)
	return m
end
LinearQuadraticModel(s::LPSolver) = LPMathProgModel(;s.options...)

function loadproblem!(m::LPMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
end

function optimize!(m::LPMathProgModel)
end

function status(m::LPMathProgModel)
end

end

