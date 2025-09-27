using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

plateDim1 = [20.0,24.0]
pointSpacing1 = 2.0

orientation1 = :up
F, V, Eb, Cb = triplate(plateDim1, pointSpacing1; orientation=orientation1, return_boundary_edges=Val(true))

# Visualization
GLMakie.closeall()

Ebs,Vbs = separate_vertices(Eb,V)
Cbs_V = simplex2vertexdata(Ebs,Cb)

fig = Figure(size=(1200,800))
ax1 = AxisGeom(fig[1, 1], title = "Triangulated mesh plate", azimuth=-pi/2, elevation=pi/2)
hp2 = meshplot!(ax1, F, V)
hp3 = edgeplot!(ax1, Ebs, Vbs; color=Cbs_V, linewidth=6, colormap=cmap)

# normalplot(ax1,F1,V1)
fig
