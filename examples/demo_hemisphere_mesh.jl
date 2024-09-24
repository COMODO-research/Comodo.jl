using Comodo
using GLMakie
using GeometryBasics
using Rotations
using Statistics
using LinearAlgebra


# Demo for building a  hemissphere from a isosaheddron .

# Example geometry for a sphere that is cut so some edges are boundary edges
n = 1# Number of refinement steps of the geodesic sphere
r = 1.0 # Sphere radius

r = 1.0
F0,V0 = hemisphere(0,r)
F1,V1 = hemisphere(1,r)
F2,V2 = hemisphere(2,r)
F3,V3 = hemisphere(3,r)

# Visualization
fig = Figure(size=(1200,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "n=0")
hp1 = poly!(ax1,GeometryBasics.Mesh(V0,F0), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "n=1")
hp2 = poly!(ax2,GeometryBasics.Mesh(V1,F1), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "n=2")
hp3 = poly!(ax3,GeometryBasics.Mesh(V2,F2), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "n=3")
hp4 = poly!(ax4,GeometryBasics.Mesh(V3,F3), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
fig