using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

#=
This is a demonstration of the capabilities of the `tridisc` function which 
generates the faces `F` and vertices `V` for a triangulated disc (circle).
=#

# Define input parameters
r = 2.0 # Radius
n = 3 # Number of refinement iterations

# Create disc mesh 
F,V = tridisc(r,n)

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A triangulated mesh of a circular domain (disc)")
# scatter!(ax1,V,markersize=10,color = :black)
hp2 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white,shading=FastShading,transparency=false)
# normalplot(ax1,F,V)
fig
