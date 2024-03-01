using Comodo
using GLMakie
using GeometryBasics

#=
This demo shows the use of `subtri` to refine triangulated meshes. Each 
original input triangle spawns 4 triangles for the regined mesh (one central 
one, and 3 at each corner). The following refinement methods are implemented: 
    
    method=:linear : This is the default method, and refines the triangles in 
    a simple linear manor through splitting. Each input edge simply obtains a 
    new mid-edge node. 
    
    method=:loop : This method features Loop-subdivision. Rather than linearly 
    splitting edges and maintaining the original coordinates, as for the linear 
    method, this method computes the new points in a special weighted sense 
    such that the surface effectively approaches a "quartic box spline". Hence 
    this method both refines and smoothes the geometry through spline 
    approximation. 
=#

## Define example input
r = 0.5 #radius
M = platonicsolid(4,r) # Get example triangulated mesh (e.g. icosahedron)
V = coordinates(M)
F = faces(M)

## Refine triangulation using `subtri` and the default :linear method

Fn1,Vn1=subtri(F,V,1) # Split once, default is same as: Fn1,Vn1=subtri(F,V,1; method="linear")
Fn2,Vn2=subtri(F,V,2) # Split twice
Fn3,Vn3=subtri(F,V,3) # Split 3 times

Fn4,Vn4=subtri(F,V,1; method=:loop) # Split once 
Fn5,Vn5=subtri(F,V,2; method=:loop) # Split twice
Fn6,Vn6=subtri(F,V,3; method=:loop) # Split 3 times

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linearly refined n=1")
wireframe!(ax1,M, linewidth=8,color=:red, overdraw=false)
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linearly refined n=2")
wireframe!(ax2,M, linewidth=8,color=:red, overdraw=false)
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linearly refined n=3")
hp1 = wireframe!(ax3,M, linewidth=8,color=:red, overdraw=false)
hp2 = poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=3,color=:white,shading=FastShading,transparency=false)
Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop-method refined n=1")
wireframe!(ax4,M, linewidth=8,color=:red, overdraw=false)
poly!(ax4,GeometryBasics.Mesh(Vn4,Fn4), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop-method refined n=2")
wireframe!(ax5,M, linewidth=8,color=:red, overdraw=false)
poly!(ax5,GeometryBasics.Mesh(Vn5,Fn5), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop-method refined n=3")
hp1 = wireframe!(ax6,M, linewidth=8,color=:red, overdraw=false)
hp2 = poly!(ax6,GeometryBasics.Mesh(Vn6,Fn6), strokewidth=3,color=:white,shading=FastShading,transparency=false)
Legend(fig[2, 4],[hp1,hp2],["Initial","Refined"])

fig
