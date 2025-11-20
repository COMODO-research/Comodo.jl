using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO

#=
This demo shows the use of `tri2quad_merge!`, which aims to create a quad 
dominated mesh by grouping adjacent triangles if the resulting quad quality is 
sufficient. It also returns any triangles which could not be grouped. 
=#

GLMakie.closeall()
for testCase = 1:3
    if testCase == 1
        w = 12.0
        h = 6.0
        pointSpacing = 1.0
        V1 = rectanglepoints(w,h,pointSpacing; dir=:acw)
        r = 1.5
        n = ceil(Int, (2*pi*r)/pointSpacing)
        V2 = circlepoints(r, n; dir=:acw)    
        VT = (V1,V2,)
        R = ([1,2],)
        P = (pointSpacing)

        VTp = deepcopy(VT) # Copy for plotting

        F,V,_ = regiontrimesh(VT,R,P; numSmoothSteps=25, gridtype=:equilateral)
    elseif testCase == 2
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)
        F, V = mergevertices(F,V)
    elseif testCase == 3
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","obj","motherChild_5k.obj")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)
        F, V = mergevertices(F,V)
    end

    Fq, Ft = tri2quad_merge!(F, V; angleThreshold=45.0, numSmoothSteps=25)

    # Visualisation
    fig = Figure(size=(1400,800))
    ax1 = AxisGeom(fig[1, 1], title="Input triangulation")
    hp1 = meshplot!(ax1, F, V, color=:white)    

    ax2 = AxisGeom(fig[1, 2], title="Quadrangulation")
    hp2 = meshplot!(ax2, Fq, V, color=:blue)    
    hp3 = meshplot!(ax2, Ft, V, color=:green)  

    Legend(fig[1, 3],[hp1, hp2, hp3],["Input triangulation", "Quads", "Triangles"])

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end