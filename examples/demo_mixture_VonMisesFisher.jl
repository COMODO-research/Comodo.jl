using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.GLMakie
using Comodo.Distributions
using Random
# Random.seed!(1) # To ensure repeatable demo behaviour

GLMakie.closeall()

μ = [[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]]
κ = [25.0, 50.0, 100.0]

spl = mixture_VonMisesFisher(μ, κ)

n = 1000
P = [Point{3,Float64}(rand(spl)) for _ in 1:n]

# Visualisation
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1])

F, V = geosphere(3, 1.0)
hp3 = meshplot!(ax1, F, V, strokewidth=0.0, color=(:white,0.5), transparency=true)
scatter!(ax1, P, markersize=10, color=:red)
screen = display(GLMakie.Screen(), fig)