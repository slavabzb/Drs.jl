#!/usr/bin/env julia

tests = [
	"transform_to_standard_form",
	"find_potential_basis",
	"solver1",
	"solver2",
	"solver3",
	"solver4"
]

for t in tests
	include("$t.jl")
end
