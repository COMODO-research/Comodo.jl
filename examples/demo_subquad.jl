using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

testCase = 2

## Define example input
if testCase == 1 # cube
    r = 1.0 #radius
    M = platonicsolid(2,r) # Get an example quadrilateral mesh (for a cube in this case)
    V = coordinates(M)
    F = faces(M)
elseif testCase == 2 # Extruded prism/cylinder with nc points
    r = 1.0
    nc = 3
    Vc = circlepoints(r,nc;dir=:cw)    
    d = norm(Vc[1]-Vc[2])
    n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))    
    F,V = extrudecurve(Vc; extent=d, direction=:positive, n=n, num_steps=2, close_loop=true,face_type=:quad)
    M = GeometryBasics.Mesh(V,F)
elseif testCase == 3 # Quadrangulated hemi-sphere
    n = 1
    r = 1.0
    F,V = quadsphere(n,r)
    VC = simplexcenter(F,V) # Finding triangle centre coordiantes
    F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces using z of central coordinates
    F,V = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed
    M = GeometryBasics.Mesh(V,F)
end
## Refine mesh using `subquad` and the default "linear" method
Fn1,Vn1=subquad(F,V,1) # Split once 
Fn2,Vn2=subquad(F,V,2) # Split twice
Fn3,Vn3=subquad(F,V,3) # Split 3 times

Fn4,Vn4=subquad(F,V,1;method=:Catmull_Clark) # Split once 
Fn5,Vn5=subquad(F,V,2;method=:Catmull_Clark) # Split twice
Fn6,Vn6=subquad(F,V,3;method=:Catmull_Clark) # Split 3 times

Fn7,Vn7=subquad(F,V,1; method=:Catmull_Clark, constrain_boundary=true) # Split once 
Fn8,Vn8=subquad(F,V,2; method=:Catmull_Clark, constrain_boundary=true) # Split twice
Fn9,Vn9=subquad(F,V,3; method=:Catmull_Clark, constrain_boundary=true) # Split 3 times


## Visualization
strokewidth1 = 1
lineWidth = 4
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linear, n=1")
wireframe!(ax1,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linearly, n=2")
wireframe!(ax2,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Linearly, n=3")
hp1 = wireframe!(ax3,M, linewidth=lineWidth,color=:red, overdraw=false)
hp2 = poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)
Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Catmull_Clark, n=1")
wireframe!(ax4,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax4,GeometryBasics.Mesh(Vn4,Fn4), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Catmull_Clark, n=2")
wireframe!(ax5,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax5,GeometryBasics.Mesh(Vn5,Fn5), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Catmull_Clark, n=3")
hp1 = wireframe!(ax6,M, linewidth=lineWidth,color=:red, overdraw=false)
hp2 = poly!(ax6,GeometryBasics.Mesh(Vn6,Fn6), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)
Legend(fig[2, 4],[hp1,hp2],["Initial","Refined"])

ax4 = Axis3(fig[3, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Catmull_Clark, constrained boundary, n=1")
wireframe!(ax4,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax4,GeometryBasics.Mesh(Vn7,Fn7), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax5 = Axis3(fig[3, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Catmull_Clark, constrained boundary, n=2")
wireframe!(ax5,M, linewidth=lineWidth,color=:red, overdraw=false)
poly!(ax5,GeometryBasics.Mesh(Vn8,Fn8), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)

ax6 = Axis3(fig[3, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Catmull_Clark, constrained boundary, n=3")
hp1 = wireframe!(ax6,M, linewidth=lineWidth,color=:red, overdraw=false)
hp2 = poly!(ax6,GeometryBasics.Mesh(Vn9,Fn9), strokewidth=strokewidth1,color=:white,shading=FastShading,transparency=false)
Legend(fig[3, 4],[hp1,hp2],["Initial","Refined"])

fig