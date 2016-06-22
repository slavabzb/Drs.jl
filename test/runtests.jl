using Base.Test

tests = ["test"]

for test in tests
	include("$(test).jl")
end

