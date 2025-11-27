using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO

#=
This demo shows the use of `tri2quad_merge_split!`. The function first uses 
`tri2quad_merge`, which the aims to create a quad dominated mesh, next the 
remaining triangles are each split into 3 quads and each quad is split in 4, 
resulting in a pure quadrangulation. 
See also `tri2quad_merge`. 
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
        fileName_mesh = joinpath(comododir(),"assets","obj","motherChild_2k.obj")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)
        F, V = mergevertices(F,V)
    end
    
    Fq, Vq, indInitial = tri2quad_merge_split!(F, V; angleThreshold=45.0, numSmoothSteps=25)
    
    # Visualisation
    fig = Figure(size=(1400,800))
    ax1 = AxisGeom(fig[1, 1], title="Input triangulation")
    hp1 = meshplot!(ax1, F, V, color=:white, strokewidth=1.0)    

    ax2 = AxisGeom(fig[1, 2], title="Quadrangulation")
    hp2 = meshplot!(ax2, Fq, Vq, color=:lightgreen, strokewidth=1.0)    
    
    Legend(fig[1, 3],[hp1, hp2],["Input triangulation", "Quads"])

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end