#!/usr/bin/env julia

DIR = dirname(@__FILE__())

for test in readdir(DIR)
	if ismatch(r"^test_.+\.jl", test)
		include("$(test)")
	end
end
