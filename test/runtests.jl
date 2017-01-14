using Base.Test

DIR = dirname(@__FILE__())

for test in readdir(DIR)
	if ismatch(r"^(?(?=runtests)$|.+\.jl$)", test)
		include("$(test)")
	end
end

