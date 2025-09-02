using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

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

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1 # icosahedron
        r = 0.5 #radius
        F,V = platonicsolid(4,r)     
    elseif testCase == 2 # Extruded prism/cylinder with nc points
        r = 1.0
        nc = 3
        Vc = circlepoints(r,nc;dir=:cw)    
        d = norm(Vc[1]-Vc[2])        
        F,V = extrudecurve(Vc; extent=d, direction=:positive, num_steps=2, close_loop=true,face_type=:backslash)    
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

    ax1 = AxisGeom(fig[1, 1], title = "Linear, n=1")
    edgeplot!(ax1, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax1, Fn1, Vn1, strokewidth=strokewidth1)

    ax2 = AxisGeom(fig[1, 2], title = "Linear, n=2")
    edgeplot!(ax2, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax2, Fn2, Vn2, strokewidth=strokewidth1)

    ax3 = AxisGeom(fig[1, 3], title = "Linear, n=3")
    hp1 = edgeplot!(ax3, F, V, linewidth=lineWidth, color=:red)
    hp2 = meshplot!(ax3, Fn3, Vn3, strokewidth=strokewidth1)
    Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

    ax4 = AxisGeom(fig[2, 1], title = "Loop, n=1")
    edgeplot!(ax4, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax4, Fn4, Vn4, strokewidth=strokewidth1)

    ax5 = AxisGeom(fig[2, 2], title = "Loop, n=2")
    edgeplot!(ax5, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax5, Fn5, Vn5, strokewidth=strokewidth1)

    ax6 = AxisGeom(fig[2, 3], title = "Loop, n=3")
    hp1 = edgeplot!(ax6, F, V, linewidth=lineWidth, color=:red)
    hp2 = meshplot!(ax6, Fn6, Vn6, strokewidth=strokewidth1)
    Legend(fig[2, 4],[hp1,hp2],["Initial","Refined"])

    ax4 = AxisGeom(fig[3, 1], title = "Loop, constrained boundary, n=1")
    edgeplot!(ax4, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax4, Fn7, Vn7, strokewidth=strokewidth1)

    ax5 = AxisGeom(fig[3, 2], title = "Loop, constrained boundary, n=2")
    edgeplot!(ax5, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax5, Fn8, Vn8, strokewidth=strokewidth1)

    ax6 = AxisGeom(fig[3, 3], title = "Loop, constrained boundary, n=3")
    hp1 = edgeplot!(ax6, F, V, linewidth=lineWidth, color=:red)
    hp2 = meshplot!(ax6, Fn9, Vn9, strokewidth=strokewidth1)
    Legend(fig[3, 4],[hp1,hp2],["Initial","Refined"])

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end