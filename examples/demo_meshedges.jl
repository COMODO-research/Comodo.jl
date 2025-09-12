using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()
testCase = 3
# for testCase = 1:2
    if testCase == 1
        F = [TriangleFace{Int64}(1,2,3)]
        V = [Point{3,Float64}(0.0, 0.0, 0.0), 
             Point{3,Float64}(1.0, 0.0, 0.0), 
             Point{3,Float64}(1.0, 1.0, 0.0)]
    elseif testCase == 2
        F = [QuadFace{Int64}(1,2,3,4)]
        V = [Point{3,Float64}(0.0, 0.0, 0.0), 
             Point{3,Float64}(1.0, 0.0, 0.0), 
             Point{3,Float64}(1.0, 1.0, 0.0),
             Point{3,Float64}(0.0, 1.0, 0.0)] 
    elseif testCase == 3
        plateDim = [5.0, 3.0]
        plateElem = [5,3]
        orientation = :up
        F, V = quadplate(plateDim, plateElem; orientation=orientation)    
    elseif testCase == 4    
        r = 1.0
        n = 0
        F,V = geosphere(n,r)
    elseif testCase == 5
        pointSpacing = 1.0
        boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        F,V,C = quadbox(boxDim,boxEl)
    end

    edgeSet = meshedges(F; unique_only=false, edgetypes=[1,2,3,4])

    Eu,indReverse = gunique(edgeSet; return_unique=Val(true), return_inverse=Val(true), sort_entries=true)
    i_e = 4
    ind1 = indReverse[(i_e-1)*length(F)+1:i_e*length(F)]

    ## Visualization
    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Mesh edges")
    hp1 = meshplot!(ax1, F, V, strokewidth=5.0)
    hp2 = edgeplot!(ax1, Eu[ind1], V, color=:red, linewidth=3)

    Legend(fig[1,2], [hp1, hp2], ["Mesh", "Edges"])
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
# end