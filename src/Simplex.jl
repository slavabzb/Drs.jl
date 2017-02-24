__precompile__()

module Simplex

export CHUZR,
  BTRAN,
  PRICE

function CHUZR(invB, basic_vars)
	r, c = size(invB)
	e = eye(r)
	rations = [basic_vars[i] / norm(invB * e[:,i]) for i in 1:r]
	P = collect(take(sortperm(rations), nprocs()))
end

function BTRAN(basic_vars, p)
  println("BTRAN pid $(myid())")
  r, c = size(basic_vars)
	e = eye(r)
	pi = e[:,p]' * basic_vars
  println("BTRAN pi $pi")
  pi
end

function PRICE(N, pi)
  println("PRICE pid $(myid())")
  pivotal_row = pi * N
  println("PRICE pivotal_row $pivotal_row")
  pivotal_row
end

end
