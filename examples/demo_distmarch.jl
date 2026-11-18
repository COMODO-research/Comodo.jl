using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

for testCase = 1:7
    if testCase == 1
        # Example geometry
        s=1.0
        V=Vector{GeometryBasics.Point{3, Float64}}(undef, 3)
        V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
        V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
        V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)

        F = [TriangleFace{Int}(1,2,3)]
        # F, V = subtri(F,V,2)

        indStart = [1]
    elseif testCase==2
        # Bowtie mesh
        s=1.0
        V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
        V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
        V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
        V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)
        V[4 ] = GeometryBasics.Point{3, Float64}( 2*s,    s, 0.0)
        V[5 ] = GeometryBasics.Point{3, Float64}( 2*s,    -s, 0.0)

        F = [TriangleFace{Int}(1,2,3), TriangleFace{Int}(5,4,3)]
        # F, V = subtri(F,V,2)

        indStart = [1]
    elseif testCase==3
        r = 1.0
        F,V = geosphere(5,r) 
        z = [v[3] for v ∈ V]
        indStart = [findmax(z)[2]]   
    elseif testCase==4
        r = 1.0
        F,V = quadsphere(r, 0.025)    
        z = [v[3] for v ∈ V]
        indStart = [findmax(z)[2]]        
    elseif testCase==5
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F, V, _ = mergevertices(F,V)
        F, V = subtri(F,V,1;method=:Loop)
        z = [v[3] for v ∈ V]
        indStart = [findmax(z)[2]]
    elseif testCase==6
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
        z = [v[3] for v ∈ V]
        indStart = [findmax(z)[2]]
    elseif testCase==7 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","femur.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
        F, V = subtri(F,V,1;method=:Loop)
        z = [v[3] for v ∈ V]
        indStart = [findmax(z)[2]]
    end

    d,dd,l = distmarch(F, V, indStart)

    ## Visualization
    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Distances")
    hp1 = meshplot!(ax1, F, V, color=d, strokewidth=0.1, colormap=Reverse(:Spectral))
    scatter!(ax1,V[indStart], color=:black, markersize=15)
    Colorbar(fig[1, 2], hp1)

    ax2 = AxisGeom(fig[1, 3], title = "Sine on distance")
    hp2 = meshplot!(ax2, F, V, color=sin.(16.0.*2.0.*pi.*(d./maximum(d))), strokewidth=0.0, colormap=:viridis)
    scatter!(ax2,V[indStart], color=:black, markersize=15)
    Colorbar(fig[1, 4], hp2)

    ax3 = AxisGeom(fig[1, 5], title = "Point regions")
    hp3 = meshplot!(ax3, F, V, color=l, strokewidth=0.0, colormap=Makie.Categorical(Reverse(:Spectral)))
    scatter!(ax3,V[indStart], color=:black, markersize=15)
    Colorbar(fig[1, 6], hp3)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end