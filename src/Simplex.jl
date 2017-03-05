__precompile__()

module Simplex

export CHUZR,
  BTRAN,
  PRICE,
  CHUZR_MI,
  CHUZC,
  UPDATE_MI,
  FTRAN1,
  FTRAN2,
  INVERT,
  UPDATE

function CHUZR(invB, basic_vars)
	r, c = size(invB)
	e = eye(r)
	rations = [basic_vars[i] / norm(invB * e[:,i]) for i in 1:r]
  println("CHUZR b/s $rations")
	collect(take(sortperm(rations, rev=true), nprocs()))
end

function BTRAN(invB, p)
  r, c = size(invB)
	e = eye(r)
	e[:,p]' * invB
end

function PRICE(N, pi)
  println("PRICE pid $(myid())")
  pi * N
end

function CHUZR_MI(P)
  splice!(P, 1)
end

function CHUZC(c, a)
  println("c $c, a $a")
  @assert length(c) == length(a)
  rations = [c[i] / a[i] for i in 1:length(c)]
  println("rations $rations")
  collect(take(sortperm(rations), 1))
end

function UPDATE_MI()
end

function FTRAN1()
end

function FTRAN2()
end

function INVERT()
end

function UPDATE()
end

end
