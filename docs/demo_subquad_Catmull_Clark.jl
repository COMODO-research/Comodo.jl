using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using LinearAlgebra

## Define example input
r = 1.0 #radius
M = platonicsolid(2,r) # Get an example quadrilateral mesh (for a cube in this case)
V = coordinates(M)
F = faces(M)

## Refine mesh using `subquad` and the Catmull-Clark method
Fn1,Vn1 = subquad(F,V,1; method = "Catmull-Clark") # Refined once 
Fn2,Vn2 = subquad(F,V,2; method = "Catmull-Clark") # Refined twice
Fn3,Vn3 = subquad(F,V,3; method = "Catmull-Clark") # Refined 3 times

Rn3 = [norm(v) for v in Vn3] 
rs3 = mean(Rn3)
d3 = Rn3.-rs3

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=1")
wireframe!(ax1,M, linewidth=8,color=:red, overdraw=false)
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=2")
wireframe!(ax2,M, linewidth=8,color=:red, overdraw=false)
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=3")
hp1 = wireframe!(ax3,M, linewidth=8,color=:red, overdraw=false)
hp2 = poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=3,color=:white,shading=FastShading,transparency=false)
Legend(fig[2, 1][1,2],[hp1,hp2],["Initial","Refined"])

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=3, distance to sphere")
hp3 = poly!(ax4,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=0.5,color=d3,shading=FastShading,transparency=false,colormap=:Spectral,limits=maximum(abs.(d3)).*(-1.0,1.0))
Colorbar(fig[2, 2][1, 2],hp3)
fig
