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
C_V1 = [v[3] for v in V1]
C_F1 = [mean(V1[f])[3] for f in F1]

function tofaceview(F, CF)
    return FaceView(CF, eltype(F).(eachindex(F)))
end    

# Visualize mesh
GLMakie.closeall()

fig = Figure(size = (1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Vertex shaded single color")
hp1 = meshplot!(ax1, F1, V1; strokewidth=1)

NF = FaceView(facenormal(F1,V1), eltype(F1).(eachindex(F1)))
M = GeometryBasics.Mesh(V1, F1, normal = NF)

ax2 = AxisGeom(fig[1, 2], title = "Face shaded single color")
hp2 = meshplot!(ax2, M; strokewidth=1)

ax3 = AxisGeom(fig[1, 3], title = "Unshaded single color")
hp3 = meshplot!(ax3, F1, V1; strokewidth=0.5, color=:green, shading=false)

ax4 = AxisGeom(fig[2, 1], title = "Transparency")
hp4 = meshplot!(ax4, F1, V1; strokewidth=2, color=(:orange,0.5), transparency=true)


ax5 = AxisGeom(fig[2, 2], title = "Vertex color data")
hp5 = meshplot!(ax5, F1, V1; strokewidth=0.5, color=C_V1)

CF = FaceView(C_F1, eltype(F1).(eachindex(F1)))
NF = FaceView(facenormal(F1,V1), eltype(F1).(eachindex(F1)))
M = GeometryBasics.Mesh(V1, F1, normal = NF, color=CF)

ax6 = AxisGeom(fig[2, 3], title = "Face color data ")
hp6 = meshplot!(ax6, M; strokewidth=0.5)

fig

