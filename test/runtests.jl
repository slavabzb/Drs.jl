#!/usr/bin/env julia

tests = [
	"transform_to_standard_form",
	"find_potential_basis"
]

for t in tests
	include("$(t).jl")
end
