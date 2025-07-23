using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

## Defining a set of points on a circle using scalar radius
r = 1.0 # Radius
n = 40 # Number of points
V1 = circlepoints(r,n; dir=:acw)

## Applying a radial function
rFun(t) = r + 0.5.*sin(3.0*t)
V2 = circlepoints(rFun,n; dir=:cw)

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A circle of " * string(n) * " points")
lines!(ax1,V1, linewidth = 5, color = :blue)
scatter!(ax1,V1,markersize=25,color = :black)
scatter!(ax1,V1[1],markersize=35,color = :red)
scatter!(ax1,V1[2],markersize=35,color = :yellow)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Radial function")
lines!(ax2,V2, linewidth = 5, color = :blue)
scatter!(ax2,V2,markersize=25,color = :black)
scatter!(ax2,V2[1],markersize=35,color = :red)
scatter!(ax2,V2[2],markersize=35,color = :yellow)

fig