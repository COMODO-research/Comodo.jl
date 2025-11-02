using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO

for testCase=1:3
    if testCase == 1
        F, V = geosphere(2, 1.0)
    elseif testCase == 2
        r = 1.0
        nf = (4,4)
        F,V = hexagonmesh(r,nf)
    elseif testCase == 3
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)
    end

    scaleFactor = 0.8
    Fl, Vl = faceedgelattice(F, V, scaleFactor)

    ## Visualisation
    GLMakie.closeall()

    markersize1 = 25
    markersize2 = 15
    linewidth1 = 2
    linewidth2 = 1

    cmap = cgrad(:viridis,3, categorical = true) 

    fig = Figure(size = (1200,800))
    ax1 = AxisGeom(fig[1,1]); 
    hp1 = meshplot!(ax1, F, V; color=:white)

    ax2 = AxisGeom(fig[1,2]);
    hp2 = meshplot!(ax2, Fl, Vl; color=:lightgreen)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end