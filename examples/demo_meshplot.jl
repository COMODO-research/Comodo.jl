using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of meshplot to visualise meshes. See also Makie's 
function poly
=#

# Example mesh
r = 1.0 # radius
n1 = 0 # Number of refinement iterations
F1,V1 = geosphere(n1,r)
c = [v[3] for v in V1]

#Visualize mesh
fig = Figure(size = (1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Visualised mesh")
hp1 = meshplot!(ax1, F1, V1; strokewidth=2)

ax2 = AxisGeom(fig[1, 2], title = "Visualised mesh")
hp2 = meshplot!(ax2, F1, V1; strokewidth=1, color=:lightgreen)

ax3 = AxisGeom(fig[2, 1], title = "Visualised mesh")
hp3 = meshplot!(ax3, F1, V1; strokewidth=0.5, color=:lightblue, stroke_depth_shift=-0.01f0)

ax3 = AxisGeom(fig[2, 2], title = "Visualised mesh")
hp3 = meshplot!(ax3, F1, V1; strokewidth=0.5, color=c, stroke_depth_shift=-0.01f0)

fig