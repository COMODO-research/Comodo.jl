using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using FileIO

c = RGBA(0.5803921568627451, 0.3411764705882353, 0.9215686274509803,0.9)

#=
This demo shows the use of the dualclag function. 
=#

GLMakie.closeall()
for testCase = 1:12
    if testCase == 1 # SIngle triangle
        n = 6
        V = circlepoints(1.0,n; dir=:acw)
        # push!(V,eltype(V)(0.0,0.0,0.0))
        F = [NgonFace{n,Int}(1:n)]
    elseif testCase == 2 # n-slice pizza triangle set
        n = 6
        V = circlepoints(1.0,n; dir=:acw)
        push!(V,eltype(V)(0.0,0.0,0.0))
        F = [TriangleFace(i,mod1(i+1,n),n+1) for i in 1:n]
    elseif testCase == 3
        F,V = tetrahedron(√3)         
    elseif testCase == 4
        F,V = cube(√3)    
    elseif testCase == 5
        F,V = icosahedron(√3)         
    elseif testCase == 6
        F,V = dodecahedron(1.0) 
    elseif testCase == 7 
        F,V = geosphere(2,1.0; method=:Loop)
    elseif testCase == 8
        F,V = subquadsphere(2,1.0)
    elseif testCase == 9 
        F,V = geosphere(2,1.0; method=:Loop)
        B = [v[3]>0.0 for v in V]
        BF = [all(B[f]) for f in F]
        F = F[BF]
        F,V = remove_unused_vertices(F,V)
    elseif testCase == 10
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V = mergevertices(F,V)
    elseif testCase == 11
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V = mergevertices(F,V)
        B = [v[3]>-10 for v in V]
        BF = [all(B[f]) for f in F]
        F = F[BF]
        F,V = remove_unused_vertices(F,V)   
    elseif testCase == 12
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","obj","motherChild_5k.obj")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M)) 
    end

    F_dual, V_dual = meshdual(F,V)

    # Visualisation
    fig = Figure(size = (1200,1200))
    ax1 = AxisGeom(fig[1, 1], title = "dual surface")
    hp1 = meshplot!(ax1, F, V; color=(:white,0.6f0), transparency=true, strokewidth=0.25f0)
    hp2 = meshplot!(ax1, GeometryBasics.Mesh(V_dual,F_dual); transparency=true, color=c, depth_shift = -0.02f0, stroke_depth_shift = -0.03f0, strokewidth=2.0f0)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end