#!/usr/bin/env julia

tests = [
	"transform_to_standard_form"
]

for t in tests
	include("$(t).jl")
end
