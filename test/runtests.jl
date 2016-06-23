using Base.Test

solvers = Any[]

if isdir(Pkg.dir("LPLib"))
	importall LPLib.LPSolverInterface
	push!(solvers, LPSolver)
end

for test in readdir("test")
	if ismatch(r"^(?(?=runtests)$|.+\.jl$)", test)
		include("$(test)")
	end
end

