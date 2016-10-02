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
println("loadproblem")
println("m: ", m)
println("A: ", A)
println("collb: ", collb)
println("obj: ", obj)
println("rowlb: ", rowlb)
println("rowub: ", rowub)
println("sense: ", sense)
end

function optimize!(m::LPMathProgModel)
println("optimize")
end

function status(m::LPMathProgModel)
println("status")
end

end

