using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra

#=
This demo shows the use of the `basisGramSchmidt` function. This function can be
used to obtain a set of mutually orthogonal basis vectors from a set of input 
vectors. 
=#

# Define non-orthogonal vector triplet 
Rz = RotXYZ(    0.0, 0.0, 0.25*π)
Rx = RotXYZ(-0.25*π,    0.0, 0.0)

u = Rz*Vec{3,Float64}(1.5, 0.0, 0.0)
v = Vec{3,Float64}(0.0, 2.0, 0.0)
w = Rx*Vec{3,Float64}(0.0, 0.0, 2.5)

A = [u,v,w]

E = basisGramSchmidt(A)

## Visualization
GLMakie.closeall()
O = zeros(Vec{3,Float64},3)
fig = Figure(size=(1800,800))
ax1 = AxisGeom(fig[1, 1])
hp1 = arrows3d!(ax1, O, A, color=[(:red, 0.25), (:green, 0.25), (:blue, 0.25)], markerscale=1.0, transparency=true)
hp2 = arrows3d!(ax1, O, E, color=[:red, :green, :blue], markerscale=1.0)

Legend(fig[1, 2],[hp1,hp2],["Vectors", "Gram-Schmidt basis"])

fig