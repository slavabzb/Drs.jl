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
	p = collect(take(sortperm(rations, rev=true), 1))
  return p[1]
end

function BTRAN(invB, p)
  r, c = size(invB)
	e = eye(r)
	return e[:,p]' * invB
end

function PRICE(N, pi)
  println("PRICE pid $(myid())")
  return pi * N
end

function CHUZR_MI(P)
  splice!(P, 1)
end

function CHUZC(c, a, nonbasis)
  cn = c[nonbasis]
  @assert length(cn) == length(a)
  rations = [cn[i] / a[i] for i in 1:length(a)]
  println("rations $rations")
  q = collect(take(sortperm(rations), 1))
  q = q[1]
  b = cn[q] / a[q]
  c[nonbasis] -= b * a'
  return q
end

function UPDATE_MI()
end

function FTRAN1(A, invB, basic_vars, pivotal_row, p, q)
  alpha = (basic_vars[p] / pivotal_row[q])[1]
  aq = invB * A[:,q]
  basic_vars[:] -= alpha * aq
end

function FTRAN2()
end

function INVERT()
end

function UPDATE()
end

end
