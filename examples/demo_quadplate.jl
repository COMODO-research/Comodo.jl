using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

plateDim = [20.0,24.0]
plateElem = [11,16]
orientation = :up
F, V, Eb, Cb = quadplate(plateDim, plateElem; orientation=orientation)

## Visualize mesh
GLMakie.closeall()

cmap = Makie.Categorical(:Spectral) 

Ebs,Vbs = separate_vertices(Eb,V)
Cbs_V = simplex2vertexdata(Ebs,Cb)

fig = Figure(size = (1200,800))
ax1 = AxisGeom(fig[1, 1], title = "Quadrilateral mesh plate")
hp2 = meshplot!(ax1, F, V)
# scatter!(ax1, V[indBoundaryNodes], markersize=15, color=:red)
hp3 = edgeplot!(ax1, Ebs, Vbs; color=Cbs_V, linewidth=6, colormap=cmap)
Colorbar(fig[1, 2], hp3)
fig