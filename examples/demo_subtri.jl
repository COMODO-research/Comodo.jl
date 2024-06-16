using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

#=
This demo shows the use of `subtri` to refine triangulated meshes. Each 
original input triangle spawns 4 triangles for the refined mesh (one central 
one, and 3 at each corner). The following refinement methods are implemented: 
    
    method=:linear : This is the default method, and refines the triangles in 
    a simple linear manor through splitting. Each input edge simply obtains a 
    new mid-edge node. 
    
    method=:Loop : This method features Loop-subdivision. Rather than linearly 
    splitting edges and maintaining the original coordinates, as for the linear 
    method, this method computes the new points in a special weighted sense 
    such that the surface effectively approaches a "quartic box spline". Hence 
    this method both refines and smoothes the geometry through spline 
    approximation. 
=#

## Define example input
testCase = 2
if testCase == 1 # icosahedron
    r = 0.5 #radius
    M = platonicsolid(4,r) 
    V = coordinates(M)
    F = faces(M)
elseif testCase == 2 # Extruded prism/cylinder with nc points
    r = 1.0
    nc = 3
    Vc = circlepoints(r,nc;dir=:cw)    
    d = norm(Vc[1]-Vc[2])        
    F,V = extrudecurve(Vc; extent=d, direction=:positive, num_steps=2, close_loop=true,face_type=:tri_slash)
    M = GeometryBasics.Mesh(V,F)
end

## Refine triangulation using `subtri` and the default :linear method

Fn1,Vn1=subtri(F,V,1) # Split once, default is same as: Fn1,Vn1=subtri(F,V,1; method="linear")
Fn2,Vn2=subtri(F,V,2) # Split twice
Fn3,Vn3=subtri(F,V,3) # Split 3 times

Fn4,Vn4=subtri(F,V,1; method=:Loop) # Split once 
Fn5,Vn5=subtri(F,V,2; method=:Loop) # Split twice
Fn6,Vn6=subtri(F,V,3; method=:Loop) # Split 3 times

Fn7,Vn7=subtri(F,V,1; method=:Loop, constrain_boundary=true) # Split once 
Fn8,Vn8=subtri(F,V,2; method=:Loop, constrain_boundary=true) # Split twice
Fn9,Vn9=subtri(F,V,3; method=:Loop, constrain_boundary=true) # Split 3 times


## Visualization
strokewidth1 = 1
lineWidth = 4
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linear, n=1")
wireframe!(ax1,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linear, n=2")
wireframe!(ax2,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linear, n=3")
hp1 = wireframe!(ax3,M, linewidth=lineWidth,color=:red, overdraw=false)
hp2 = poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)
Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop, n=1")
wireframe!(ax4,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax4,GeometryBasics.Mesh(Vn4,Fn4), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop, n=2")
wireframe!(ax5,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax5,GeometryBasics.Mesh(Vn5,Fn5), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop, n=3")
hp1 = wireframe!(ax6,M, linewidth=lineWidth,color=:red, overdraw=false)
hp2 = poly!(ax6,GeometryBasics.Mesh(Vn6,Fn6), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)
Legend(fig[2, 4],[hp1,hp2],["Initial","Refined"])

ax4 = Axis3(fig[3, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop, constrained boundary, n=1")
wireframe!(ax4,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax4,GeometryBasics.Mesh(Vn7,Fn7), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax5 = Axis3(fig[3, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop, constrained boundary, n=2")
wireframe!(ax5,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax5,GeometryBasics.Mesh(Vn8,Fn8), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax6 = Axis3(fig[3, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Loop, constrained boundary, n=3")
hp1 = wireframe!(ax6,M, linewidth=lineWidth,color=:red, overdraw=false)
hp2 = poly!(ax6,GeometryBasics.Mesh(Vn9,Fn9), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)
Legend(fig[3, 4],[hp1,hp2],["Initial","Refined"])

fig
