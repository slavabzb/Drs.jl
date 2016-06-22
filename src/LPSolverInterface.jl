module LPSolverInterface

importall MathProgBase.SolverInterface

export LPSolver

type LPMathProgModel <: AbstractLinearQuadraticModel
end

immutable LPSolver <: AbstractMathProgSolver
    options
end
LPSolver(;kwargs...) = LPSolver(kwargs)

end
