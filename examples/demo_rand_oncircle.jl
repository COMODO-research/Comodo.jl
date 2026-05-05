using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.GLMakie
using Comodo.Distributions
using Comodo.Rotations 
using Random

GLMakie.closeall()

r = 1.0 # Sphere radius
N = 100 # Number of points 
P = rand_oncircle(r, N) # Points inside the sphere

# Visualisation
fig = Figure(size=(800,800))
ax1 = Axis(fig[1, 1], aspect = DataAspect(), title="random points on circle")

Vc = circlepoints(r, 1000) 
lines!(ax1, Vc, linewidth=3, color=:red)
scatter!(ax1, P, markersize=15, color=:black)
screen = display(GLMakie.Screen(), fig)