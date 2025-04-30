using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using FileIO

using Comodo.Statistics
using Comodo.LinearAlgebra

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1    
        fileName_mesh = joinpath(comododir(),"assets","stl","femur.stl")
        M = load(fileName_mesh)
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = [Point{3,Float64}(p) for p in coordinates(M)]
        F,V  = mergevertices(F,V)
        M = GeometryBasics.Mesh(V,F)
        indFace = 1
    elseif testCase == 2
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")    
        M = load(fileName_mesh)
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = [Point{3,Float64}(p) for p in coordinates(M)]
        F,V  = mergevertices(F,V)
        M = GeometryBasics.Mesh(V,F)
        indFace = 1
    elseif testCase == 3
        n = 4 
        r = 1.0
        F,V = geosphere(4,2.0)
        M = GeometryBasics.Mesh(V,F)
        indFace = 1
    end

    P_on1 = faceinteriorpoint(F,V, indFace; w=0.0)
    P_in = faceinteriorpoint(F,V, indFace; w=0.5)
    P_on2 = faceinteriorpoint(F,V, indFace; w=1.0)

    ## Visualization
    c = RGBA(0.5, 0.5, 0.5,0.25)

    markerSize = 15

    fig = Figure(size=(800,800))
    ax1 = LScene(fig[1,1]) # ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :ray, triSide=1""")
    hp1 = poly!(ax1,M,color=c, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.25)

    hp2 = scatter!(ax1,P_on1,color=:red,markersize=markerSize)
    hp3 =scatter!(ax1,P_in,color=:green,markersize=markerSize)
    hp4 =scatter!(ax1,P_on2,color=:blue,markersize=markerSize)
    hp5 =lines!(ax1,[P_on1, P_on2],color=:yellow)
    Legend(fig[1, 2],[hp2,hp3,hp4],["on1", "in", "on2"])

    stepRange = 1:1:length(F)
    hSlider = Slider(fig[2, :], range = stepRange, startvalue = indFace,linewidth=30)

    on(hSlider.value) do sliderVal 
        P_on1 = faceinteriorpoint(F,V, sliderVal; w=0.0)
        P_in = faceinteriorpoint(F,V, sliderVal; w=0.5)
        P_on2 = faceinteriorpoint(F,V, sliderVal; w=1.0)
        hp2[1] = P_on1
        hp3[1] = P_in    
        hp4[1] = P_on2
        hp5[1] = [P_on1,P_on2]
    end

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end