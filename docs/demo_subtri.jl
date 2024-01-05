using Gibbon
using GLMakie
using GeometryBasics

r = 0.5 #radius
n = 3

M = platonicsolid(4,r)
V = coordinates(M)
F = faces(M)

Fn,Vn = subtri(F,V,n)


# Visualize mesh
fig = Figure(size=(800,800))

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined Icosahedron")

hp1=wireframe!(ax,M,linewidth=8,color=:red, overdraw=false)
hp2=poly!(ax,GeometryBasics.Mesh(Vn,Fn),strokewidth=1,color=:white, shading = FastShading)

Legend(fig[1, 2],[hp1,hp2],["Initial","Refined"])
fig