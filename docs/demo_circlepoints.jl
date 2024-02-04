using Gibbon
using GLMakie
using GeometryBasics
using LinearAlgebra

r = 1.0 # Radius
n = 25 # Number of points

V = circlepoints(r,n)

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A circle of " * string(n) * " points")
lines!(ax1,V, linewidth = 5, color = :blue)
scatter!(V,markersize=25,color = :black)
fig