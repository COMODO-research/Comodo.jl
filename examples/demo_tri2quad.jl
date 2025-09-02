using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

#=
This demo is for the `tri2quad` function which can be used to convert triangles 
to quadrilaterals. Quads are created at each triangle vertex and are formed by 
connecting the original vertices to new mid-edge and mid-face vertices.  
=#

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1 # Single triangle
        V = [Point{3,Float64}(0.0,0.0,0.0),Point{3,Float64}(1.0,0.0,0.0),Point{3,Float64}(1.0,1.0,0.0)]
        F = [TriangleFace{Int}(1,2,3)]
    elseif testCase == 2 # Single triangle refined once, so 4 triangles
        V = [Point{3,Float64}(0.0,0.0,0.0),Point{3,Float64}(1.0,0.0,0.0),Point{3,Float64}(0.0,1.0,0.0),Point{3,Float64}(1.0,1.0,0.0)]
        F = [TriangleFace{Int}(1,2,3),TriangleFace{Int}(2,4,3)]
        F,V = subtri(F,V,1)
    elseif testCase == 3 # tetrahedron
        r = 0.5 #radius
        F,V = platonicsolid(1,r)     
    elseif testCase == 4 # geosphere
        r = 0.5 #radius
        F,V = geosphere(1,r) 
    elseif testCase == 5 # bunny
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)
        
        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        
        F,V = mergevertices(F,V)
    end


    Fq,Vq = tri2quad(F,V; method=:split)
    Fr,Vr = tri2quad(F,V; method=:rhombic)

    ## Visualization

    Fp,Vp = separate_vertices(F,V)
    Fqp,Vqp = separate_vertices(Fq,Vq)
    Frp,Vrp = separate_vertices(Fr,Vr)

    strokewidth1 = 1
    lineWidth = 4
    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Input triangulation")
    meshplot!(ax1, Fp, Vp, strokewidth=strokewidth1)

    ax2 = AxisGeom(fig[1, 2], title = "split method")    
    hp2 = meshplot!(ax2, Fqp, Vqp, strokewidth=strokewidth1)
    # scatter!(ax2,Vq,markersize=25,color=:green)
    # normalplot(ax2,Fq,Vq)

    ax3 = AxisGeom(fig[1, 3], title = "rhombic method")    
    hp2 = meshplot!(ax3, Frp, Vrp, strokewidth=strokewidth1)
    # scatter!(ax3,Vr,markersize=25,color=:green)
    # normalplot(ax3,Fr,Vr)
    # Legend(fig[1, 3],[hp1,hp2],["Input triangulation","Output quadrangulation"])
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end