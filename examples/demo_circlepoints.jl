using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

## Defining a set of points on a circle using scalar radius
r = 1.0 # Radius
n = 40 # Number of points
V1 = circlepoints(r,n)

## Applying a radial function
rFun(t) = r + 0.5.*sin(3*t)
V2 = circlepoints(rFun,n)

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A circle of " * string(n) * " points")
lines!(ax1,V1, linewidth = 5, color = :blue)
scatter!(ax1,V1,markersize=25,color = :black)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Radial function")
lines!(ax2,V2, linewidth = 5, color = :blue)
scatter!(ax2,V2,markersize=25,color = :black)

fig