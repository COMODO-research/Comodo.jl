using Gibbon
using GLMakie
using GeometryBasics

## Define example input
r = 0.5 #radius
M = platonicsolid(4,r) # Get example triangulated mesh (e.g. icosahedron)
V = coordinates(M)
F = faces(M)

## Refine triangulation using `subTri`
n = 3 # Number of refinement steps
Fn,Vn = subTri(F,V,n) # Subdivide triangulation linearly 

## Visualization
fig = Figure(size=(1600,800))
ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined Icosahedron")

hp1 = wireframe!(ax,M,linewidth=8,color=:red, overdraw=false)
hp2 = poly!(ax,GeometryBasics.Mesh(Vn,Fn),strokewidth=1,color=:white, shading = FastShading)
Legend(fig[1, 2],[hp1,hp2],["Initial","Refined"])

fig