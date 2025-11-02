using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.LinearAlgebra

## Define example input
r = 1.0 #radius
F,V = platonicsolid(2,r) # Get an example quadrilateral mesh (for a cube in this case)

## Refine mesh using `subquad` and the Catmull-Clark method
Fn1,Vn1 = subquad(F,V,1; method = :Catmull_Clark) # Refined once 
Fn2,Vn2 = subquad(F,V,2; method = :Catmull_Clark) # Refined twice
Fn3,Vn3 = subquad(F,V,3; method = :Catmull_Clark) # Refined 3 times

Rn3 = [norm(v) for v in Vn3] 
rs3 = mean(Rn3)
d3 = Rn3.-rs3

## Visualization
GLMakie.closeall()

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Refined n=1")
edgeplot!(ax1, F, V, linewidth=8,color=:red)
meshplot!(ax1, Fn1, Vn1, strokewidth=3)

ax2 = AxisGeom(fig[1, 2], title = "Refined n=2")
edgeplot!(ax2, F, V, linewidth=8,color=:red)
meshplot!(ax2, Fn2, Vn2, strokewidth=3)

ax3 = AxisGeom(fig[2, 1], title = "Refined n=3")
hp1 = edgeplot!(ax3, F, V, linewidth=8,color=:red)
hp2 = meshplot!(ax3, Fn3, Vn3, strokewidth=3)
Legend(fig[2, 1][1,2],[hp1,hp2],["Initial","Refined"])

ax4 = AxisGeom(fig[2, 2], title = "Refined n=3, distance to sphere")
hp3 = meshplot!(ax4, Fn3, Vn3, strokewidth=0.5, color=d3, colormap=:Spectral)
Colorbar(fig[2, 2][1, 2],hp3)
fig