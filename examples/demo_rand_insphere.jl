using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.GLMakie
using Comodo.Distributions
using Comodo.Rotations 
using Random

GLMakie.closeall()

r = 1.0 # Sphere radius
N = 1500 # Number of points 
P = rand_insphere(r, N) # Points inside the sphere

# Visualisation
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title="random points in sphere")

F, V = geosphere(3, r)
hp3 = meshplot!(ax1, F, V, strokewidth=0.0, color=(:white,0.5), transparency=true)
scatter!(ax1, P, markersize=5, color=:black)
screen = display(GLMakie.Screen(), fig)