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
GLMakie.closeall()

fig = Figure(size=(1200,800))

ax1 = AxisGeom(fig[1, 1], title = "face_type=:tri, n=0")
hp1 = meshplot!(ax1, F1, V1)

ax2 = AxisGeom(fig[1, 2], title = "face_type=:tri, n=1")
hp2 = meshplot!(ax2, F2, V2)
ax3 = AxisGeom(fig[1, 3], title = "face_type=:tri, n=2")
hp3 = meshplot!(ax3, F3, V3)
ax4 = AxisGeom(fig[1, 4], title = "face_type=:tri, n=3")
hp4 = meshplot!(ax4, F4, V4)

ax5 = AxisGeom(fig[2, 1], title = "face_type=:quad, n=0")
hp5 = meshplot!(ax5, F5, V5)
# normalplot(ax5,F5,V5)
ax6 = AxisGeom(fig[2, 2], title = "face_type=:quad, n=1")
hp6 = meshplot!(ax6, F6, V6)
ax7 = AxisGeom(fig[2, 3], title = "face_type=:quad, n=2")
hp7 = meshplot!(ax7, F7, V7)

ax8 = AxisGeom(fig[2, 4], title = "face_type=:quad, n=3")
F8s,V8s = separate_vertices(F8,V8) # Give each face its own point set 
C8s = simplex2vertexdata(F8s,C8) # Convert face color data to vertex color data 
hp8 = meshplot!(ax8, F8s, V8s; color=C8s)

fig