using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `gridpoints_equilateral` to generate an equilateralt
triangle based point distribution. 
=#

pointSpacing = 1
xSpan = [-5,5]
ySpan = [-5,5]

V1 = gridpoints_equilateral(xSpan,ySpan,pointSpacing)
F2,V2 = gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = Val(true))
F3,V3 = gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = Val(true), rectangular = Val(true))

# Visualisation
fig = Figure(size=(1600,400))

ax1 = Axis3(fig[1, 1],aspect = :data,title="An equilateral triangle based point distribution",azimuth=-pi/2,elevation=pi/2)
hp1 = scatter!(ax1, V1,markersize=12,color=:black)

ax2 = Axis3(fig[1, 2],aspect = :data,title="An equilateral triangle mesh",azimuth=-pi/2,elevation=pi/2)
# hp1 = scatter!(ax2, V2,markersize=12,color=:black)
hp2 = poly!(ax2,GeometryBasics.Mesh(V2,F2), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax3 = Axis3(fig[1, 3],aspect = :data,title="Forced to rectangular domain",azimuth=-pi/2,elevation=pi/2)
# hp1 = scatter!(ax3, V3,markersize=12,color=:black)
hp2 = poly!(ax3,GeometryBasics.Mesh(V3,F3), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

fig
