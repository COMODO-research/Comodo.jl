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
F4, V4, Eb, Cb = gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = Val(true), rectangular = Val(true), return_boundary_edges=Val(true))

# Visualisation
GLMakie.closeall()

fig = Figure(size=(1600,400))

ax1 = AxisGeom(fig[1, 1], title="An equilateral triangle based point distribution", azimuth=-pi/2,elevation=pi/2)
hp1 = scatter!(ax1, V1, markersize=12,color=:black)

ax2 = AxisGeom(fig[1, 2], title="An equilateral triangle mesh", azimuth=-pi/2,elevation=pi/2)
hp1 = scatter!(ax2, V2, markersize=12,color=:black)
hp2 = meshplot!(ax2, F2, V2)

ax3 = AxisGeom(fig[1, 3], title="Forced to rectangular domain", azimuth=-pi/2,elevation=pi/2)
hp1 = scatter!(ax3, V3, markersize=12,color=:black)
hp2 = meshplot!(ax3, F3, V3)

fig
