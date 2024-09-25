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
F1,V1 = hemisphere(1,r; face_type=:tri)
F2,V2 = hemisphere(2,r; face_type=:tri)
F3,V3 = hemisphere(3,r; face_type=:tri)

F4,V4 = hemisphere(1,r; face_type=:quad)
F5,V5 = hemisphere(2,r; face_type=:quad)
F6,V6 = hemisphere(3,r; face_type=:quad)

# Visualization
fig = Figure(size=(1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=0")
hp1 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=2")
hp2 = poly!(ax2,GeometryBasics.Mesh(V2,F2), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=3")
hp4 = poly!(ax3,GeometryBasics.Mesh(V3,F3), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=0")
hp4 = poly!(ax4,GeometryBasics.Mesh(V4,F4), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=2")
hp5 = poly!(ax5,GeometryBasics.Mesh(V5,F5), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=3")
hp6 = poly!(ax6,GeometryBasics.Mesh(V6,F6), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

fig