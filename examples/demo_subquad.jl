using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

for testCase = 1:3
    ## Define example input
    if testCase == 1 # cube
        r = 1.0 #radius
        F,V = platonicsolid(2,r) # Get an example quadrilateral mesh (for a cube in this case)
        M = GeometryBasics.Mesh(V,F)
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
        VC = simplexcenter(F,V) # Finding triangle centre coordinates
        F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces using z of central coordinates
        F,V = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed
        M = GeometryBasics.Mesh(V,F)
    end
    ## Refine mesh using `subquad` and the default "linear" method
    Fn1,Vn1 = subquad(F, V, 1) # Split once 
    Fn2,Vn2 = subquad(F, V, 2) # Split twice
    Fn3,Vn3 = subquad(F, V, 3) # Split 3 times

    Fn4,Vn4 = subquad(F, V, 1; method=:Catmull_Clark) # Split once 
    Fn5,Vn5 = subquad(F, V, 2; method=:Catmull_Clark) # Split twice
    Fn6,Vn6 = subquad(F, V, 3; method=:Catmull_Clark) # Split 3 times

    Fn7,Vn7 = subquad(F, V, 1; method=:Catmull_Clark, constrain_boundary=true) # Split once 
    Fn8,Vn8 = subquad(F, V, 2; method=:Catmull_Clark, constrain_boundary=true) # Split twice
    Fn9,Vn9 = subquad(F, V, 3; method=:Catmull_Clark, constrain_boundary=true) # Split 3 times

    ## Visualization    
    lineWidth = 4
    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1], title = "Linear, n=1")
    edgeplot!(ax1, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax1, Fn1, Vn1)

    ax2 = AxisGeom(fig[1, 2], title = "Linear, n=2")
    edgeplot!(ax2, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax2, Fn2, Vn2)

    ax3 = AxisGeom(fig[1, 3], title = "Linear, n=3")
    hp1 = edgeplot!(ax3, F, V, linewidth=lineWidth, color=:red)
    hp2 = meshplot!(ax3, Fn3, Vn3)
    Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

    ax4 = AxisGeom(fig[2, 1], title = "Catmull_Clark, n=1")
    edgeplot!(ax4, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax4, Fn4, Vn4)

    ax5 = AxisGeom(fig[2, 2], title = "Catmull_Clark, n=2")
    edgeplot!(ax5, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax5, Fn5, Vn5)

    ax6 = AxisGeom(fig[2, 3], title = "Catmull_Clark, n=3")
    hp1 = edgeplot!(ax6, F, V, linewidth=lineWidth, color=:red)
    hp2 = meshplot!(ax6, Fn6, Vn6)
    Legend(fig[2, 4],[hp1,hp2],["Initial","Refined"])

    ax4 = AxisGeom(fig[3, 1], title = "Catmull_Clark, constrained boundary, n=1")
    edgeplot!(ax4, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax4, Fn7, Vn7)

    ax5 = AxisGeom(fig[3, 2], title = "Catmull_Clark, constrained boundary, n=2")
    edgeplot!(ax5, F, V, linewidth=lineWidth, color=:red)
    meshplot!(ax5, Fn8, Vn8)

    ax6 = AxisGeom(fig[3, 3], title = "Catmull_Clark, constrained boundary, n=3")
    hp1 = edgeplot!(ax6, F, V, linewidth=lineWidth, color=:red)
    hp2 = meshplot!(ax6, Fn9, Vn9)
    Legend(fig[3, 4],[hp1,hp2],["Initial","Refined"])

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end