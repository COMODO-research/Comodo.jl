using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

for testCase = 1:9
    if testCase == 1
        s=1.0
        V=Vector{Point{3, Float64}}(undef,5)
        V[1 ] = Point{3, Float64}( 0.0,    s, 0.0)
        V[2 ] = Point{3, Float64}( 0.0,   -s, 0.0)
        V[3 ] = Point{3, Float64}(   s,  0.0, 0.0)

        F = Vector{TriangleFace{Int}}(undef,1)
        F[1 ] = TriangleFace{Int}(1,2,3)
        # F,V=subtri(F,V,2)
    elseif testCase==2
        s=1.0
        V=Vector{Point{3, Float64}}(undef,5)
        V[1 ] = Point{3, Float64}( 0.0,    s, 0.0)
        V[2 ] = Point{3, Float64}( 0.0,   -s, 0.0)
        V[3 ] = Point{3, Float64}(   s,  0.0, 0.0)
        V[4 ] = Point{3, Float64}( 2*s,    s, 0.0)
        V[5 ] = Point{3, Float64}( 2*s,    -s, 0.0)

        F = Vector{TriangleFace{Int}}(undef,2)
        F[1 ] = TriangleFace{Int}(1,2,3)
        F[2 ] = TriangleFace{Int}(3,4,5)
        # F,V=subtri(F,V,2)
    elseif testCase==3
        r = 1.0
        F,V = geosphere(2,r)
    elseif testCase==4
        r = 1.0
        F,V = geosphere(2,r)
        for i = 0:1:2
            for j = 0:1:2
                for k = 0:1:2
                    if i+j+k>0
                        F2,V2 = geosphere(rand(0:2,1)[1],r)
                        F2 = map(f-> f.+length(V),F2)
                        V2 = map(v-> Point{3, Float64}(i*3*r+v[1],j*3*r+v[2],k*3*r+v[3]),V2)
                        append!(F,F2)
                        append!(V,V2)
                    end
                end
            end
        end
    elseif testCase==5
        r = 1.0
        F,V = subquadsphere(2,r)
        for i = 0:1:2
            for j = 0:1:2
                for k = 0:1:2
                    if i+j+k>0
                        F2,V2 = subquadsphere(rand(0:2,1)[1],r)
                        F2 = map(f-> f.+length(V),F2)
                        V2 = map(v-> Point{3, Float64}(i*3*r+v[1],j*3*r+v[2],k*3*r+v[3]),V2)
                        append!(F,F2)
                        append!(V,V2)
                    end
                end
            end
        end
    elseif testCase==6 # Unmerged STL, each triangle is separate group
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
    elseif testCase==7 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
    elseif testCase==8 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)

        F2 = map(f-> f.+length(V),deepcopy(F))
        V2 = map(v-> Point{3, Float64}(v[1],150+v[2],v[3]),deepcopy(V))
        append!(F,F2)
        append!(V,V2)
    elseif testCase==9
        fileName_mesh = joinpath(comododir(),"assets","obj","lego_figure.obj")    
        M = load(fileName_mesh)
        F = tofaces(faces(M))
        V = topoints(coordinates(M))    
    end

    C = meshgroup(F)
    numGroups = maximum(C)

    # c = cgrad(:Spectral,numGroups,categorical = true)
    c = Makie.Categorical(Makie.Reverse(:Spectral))

    Fn,Vn = separate_vertices(F,V)
    Cn = simplex2vertexdata(Fn,C,Vn)

    # Visualization
    fig = Figure(size=(1200,1200))

    ax1 = AxisGeom(fig[1, 1], title = "A multi-object mesh")
    hp1 = meshplot!(ax1, F, V)

    ax2 = AxisGeom(fig[1, 2], title = "Grouping")
    hp2 = meshplot!(ax2, Fn, Vn, color=Cn, colormap=c)

    Legend(fig[1, 3],[hp1,hp2],["Mesh object","Grouped mesh"])
    Colorbar(fig[1, 4],hp2, label = "Group labelling")

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end