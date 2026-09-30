using Comodo
using Comodo.GLMakie

#=
This demo shows the use of `tetrahedron` to generate the faces and vertices of a cube. 
=#

r = 1.0
F, V = tetrahedron(r)
C = collect(1:length(F))

# Visualisation
cmap = Makie.Categorical(:Spectral) 

Fbs,Vbs = separate_vertices(F,V)
Cbs_V = simplex2vertexdata(Fbs, C)

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Tetrahedron mesh")
hp2 = meshplot!(ax1, Fbs, Vbs; strokewidth=3, color=Cbs_V, colormap=cmap)
hp3 = normalplot(ax1, Fbs, Vbs; type_flag=:face, color=:black,linewidth=3)

Colorbar(fig[1, 2], hp2)

screen = display(GLMakie.Screen(), fig)