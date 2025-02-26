using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# Demo for building a  hemissphere from an isosaheddron.

# Example geometry for a sphere that is cut so some edges are boundary edges
n = 1# Number of refinement steps of the geodesic sphere
r = 1.0 # hemisphere radius

F1,V1,C1 = hemisphere(0,r; face_type=:tri,closed=true)
F2,V2,C2 = hemisphere(1,r; face_type=:tri)
F3,V3,C3 = hemisphere(2,r; face_type=:tri,closed=true)
F4,V4,C4 = hemisphere(3,r; face_type=:tri)

F5,V5,C5 = hemisphere(0,r; face_type=:quad,closed=true)
F6,V6,C6 = hemisphere(1,r; face_type=:quad)
F7,V7,C7 = hemisphere(2,r; face_type=:quad,closed=true)
F8,V8,C8 = hemisphere(2,r; face_type=:quad,closed=true)

# Visualization
fig = Figure(size=(1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=0")
hp1 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
# normalplot(ax1,F1,V1)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=1")
hp2 = poly!(ax2,GeometryBasics.Mesh(V2,F2), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=2")
hp3 = poly!(ax3,GeometryBasics.Mesh(V3,F3), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax4 = Axis3(fig[1, 4], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:tri, n=3")
hp4 = poly!(ax4,GeometryBasics.Mesh(V4,F4), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax5 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=0")
hp5 = poly!(ax5,GeometryBasics.Mesh(V5,F5), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
# normalplot(ax5,F5,V5)
ax6 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=1")
hp6 = poly!(ax6,GeometryBasics.Mesh(V6,F6), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
ax7 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=2")
hp7 = poly!(ax7,GeometryBasics.Mesh(V7,F7), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax8 = Axis3(fig[2, 4], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face_type=:quad, n=3")
F8s,V8s = separate_vertices(F8,V8) # Give each face its own point set 
C8s = simplex2vertexdata(F8s,C8) # Convert face color data to vertex color data 
hp8 = poly!(ax8,GeometryBasics.Mesh(V8s,F8s), strokewidth=1,color=C8s, strokecolor=:black, shading = FastShading, transparency=false)

fig