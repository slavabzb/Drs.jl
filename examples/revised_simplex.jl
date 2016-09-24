# https://en.wikipedia.org/wiki/Revised_simplex_method

# Consider a linear program where

c = Float64[-2 -3 -4 0 0]'
A = Float64[3 2 1 1 0; 2 5 3 0 1]
b = Float64[10 15]'

# Let initially

B = Float64[4 5]
N = Float64[1 2 3]

# It corresponds to feasible vertex

x = Float64[0 0 0 10 15]

# At this moment

lambda = Float64[0 0]'
sn = Float64[-2 -3 -4]'

# Choose m < q <= n such that sq < 0 as the entering index

q = 3

# Then

d = Float64[1 3]'

# After the pivot operation

B = Float64[3 4]
N = Float64[1 2 5]

# Correspondingly

x = Float64[0 0 5 5 0]'
lambda = Float64[0 -4/3]'
sn = Float64[2/3 11/3 4]'

# A positive sn indicates that x is now optimal

