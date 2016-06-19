using Base.Test

tests = ["test"]

for test in tests
	@printf("Running %s ", test)
	include("$(test).jl")
	@printf("%20s", "... OK\n")
end

