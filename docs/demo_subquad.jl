using Comodo
using GLMakie
using GeometryBasics

## Define example input
r = 0.5 #radius
M = platonicsolid(2,r) # Get an example quadrilateral mesh (for a cube in this case)
V = coordinates(M)
F = faces(M)

## Refine mesh using `subquad` and the default "linear" method
Fn1,Vn1=subquad(F,V,1) # Split once 
Fn2,Vn2=subquad(F,V,2) # Split twice
Fn3,Vn3=subquad(F,V,3) # Split 3 times

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=1")
wireframe!(ax1,M, linewidth=8,color=:red, overdraw=false)
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=2")
wireframe!(ax2,M, linewidth=8,color=:red, overdraw=false)
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=3")
hp1 = wireframe!(ax3,M, linewidth=8,color=:red, overdraw=false)
hp2 = poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=3,color=:white,shading=FastShading,transparency=false)

Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

fig