using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        V = [Point{3, Float64}(-1.0, -0.5, 0.0), 
            Point{3, Float64}( 1.0, -0.5, 0.0),
            Point{3, Float64}( 0.0,  1.0, 0.0),
            Point{3, Float64}( 0.0,  0.0, 0.0),
            Point{3, Float64}( 2.0,  0.0, 0.0),
            ]

        F = [TriangleFace{Int}(1,2,4),
            TriangleFace{Int}(2,3,4),
            TriangleFace{Int}(3,1,4),
            TriangleFace{Int}(3,2,5),
            ]
    elseif testCase == 2
        # Create disc mesh 
        r = 1.0 # Radius
        nd = 2
        F, V = tridisc(r, nd)
        indReplace = 1:5:length(F)
        Fn, Vn = subtri_centre(F[indReplace], V, 1) # Make three connected 
        deleteat!(F,indReplace)
        append!(F, [TriangleFace{Int}(f.+length(V)) for f in Fn])
        append!(V,Vn)
        F, V = mergevertices(F, V)
        F, V = remove_unused_vertices(F,V)
    elseif testCase == 3
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)
        F, V = mergevertices(F, V)
        F, V = subtri_centre(F, V, 2) # Make three connected 
    end

    Fn = deepcopy(F)
    Vn = deepcopy(V)
    removethreeconnected!(Fn, Vn)

    # Visualisation
    fig = Figure(size = (800,800))
    ax1 = AxisGeom(fig[1, 1], title="Triangulation featuring 3 connected points")
    hp1 = meshplot!(ax1, F,  V, color=:white)

    ax2 = AxisGeom(fig[1, 2], title="Triangulation after interative 3 connected point removal")
    hp2 = meshplot!(ax2, Fn,  Vn, color=:white)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end