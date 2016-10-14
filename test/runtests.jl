using Base.Test

for test in readdir("test")
	if ismatch(r"^(?(?=runtests)$|.+\.jl$)", test)
		include("$(test)")
	end
end

