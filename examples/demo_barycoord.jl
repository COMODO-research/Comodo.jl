using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics

#=
This demo shows the use of `barycoord` to compute the i-th Barycentric 
coordinate `λ` for the point `p` and the input triangle `f` and vertices `V`. 
=#

# Create example triangle
f, V = equilateraltriangle(2.0) # Single equilateral triangle domain

# Compute Cartesian coordinates
pm = mean(V)

# Bary centric coordinates of central point should produce 1/3
λ1m = barycoord(V, mean(V), 1)
λ2m = barycoord(V, mean(V), 2)
λ3m = barycoord(V, mean(V), 3)

# i-th Bary centric coordinates of i-th vertices should produce 1
λ1 = barycoord(V, V[1], 1)
λ2 = barycoord(V, V[2], 2)
λ3 = barycoord(V, V[3], 3)

# i-th Bary centric coordinates of not i-th vertices should produce 0
λ1 = barycoord(V, V[1], 2)
λ2 = barycoord(V, V[2], 3)
λ3 = barycoord(V, V[3], 1)