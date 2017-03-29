#!/usr/bin/env julia

A = Float64[1 -2 -1 0 0; 0 3 4 1 0; 0 6 1 0 1]
b = Float64[0, 6, 3]
c = Float64[]

println("A $A\nb $b\nc $c")

basis = [1, 4, 5]
nonbasis = [2, 3]

B = A[:,basis]
N = A[:,nonbasis]

invB = inv(B)
println("invB $invB")

terminate = false
iter = 0
maxiter = 3

while !terminate
	println("iter $iter")
	println("basis $basis\nnonbasis $nonbasis")
	
	B = A[:,basis]
	N = A[:,nonbasis]
	println("B $B\nN $N")
	
	delta = invB[1,:]' * N
	println("delta $delta")
	
	incoming = indmin(delta)
	terminate = delta[incoming] > 0
	if terminate
		continue
	end
	println("incoming $incoming")
	
	Xk = invB * N[:,incoming]
	println("Xk $Xk")
	
	theta = b[2:end] ./ Xk[2:end]
	println("theta $theta")
	
	outgoing = indmin(theta) + 1
	println("outgoing $outgoing")
	
	rect = [invB[:,2:end] b Xk]
	println("rect $rect")
	
	# gaussian elimination
	m, n = size(rect)
	
	rect[outgoing,:] /= rect[outgoing,end]
	
	for i in 1:m
		if i != outgoing
			rect[i,:] -= rect[outgoing,:] * rect[i,end]
		end
	end
	
	println("rect $rect")
	
	b = rect[:,end-1]
	println("b $b")
	
	invB[:,2:end] = rect[:,1:end-2]
	println("invB $invB")
	
	basis[outgoing], nonbasis[incoming] = nonbasis[incoming], basis[outgoing]
	
	iter += 1
	if iter > maxiter
		terminate = true
	end
end

objval = b[basis[1]]
optimal = b[basis[2:end]]

println("objval $objval\noptimal $optimal")

