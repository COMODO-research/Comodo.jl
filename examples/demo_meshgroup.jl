using Comodo
using GLMakie
using GeometryBasics

# Example geometry
testCase = 3
if testCase == 1
    s=1.0
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)

    F = Vector{TriangleFace{Int64}}(undef,1)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    # F,V=subtri(F,V,2)
elseif testCase==2
    s=1.0
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}( 2*s,    s, 0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}( 2*s,    -s, 0.0)

    F = Vector{TriangleFace{Int64}}(undef,2)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    F[2 ] = TriangleFace{Int64}(3,4,5)
    # F,V=subtri(F,V,2)
elseif testCase==3
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
elseif testCase==4
    r = 1.0
    F,V = quadsphere(2,r)
    for i = 0:1:2
        for j = 0:1:2
            for k = 0:1:2
                if i+j+k>0
                    F2,V2 = quadsphere(rand(0:2,1)[1],r)
                    F2 = map(f-> f.+length(V),F2)
                    V2 = map(v-> Point{3, Float64}(i*3*r+v[1],j*3*r+v[2],k*3*r+v[3]),V2)
                    append!(F,F2)
                    append!(V,V2)
                end
            end
        end
    end
elseif testCase==5 # Unmerged STL, each triangle is seperate group
    # Loading a mesh
    using FileIO
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
elseif testCase==6 # Merged STL for single object
    # Loading a mesh
    using FileIO
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase==7 # Merged STL for single object
    # Loading a mesh
    using FileIO
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
end

C = meshgroup(F)
numGroups = maximum(C)

c = cgrad(:Spectral,numGroups,categorical = true)

Fn,Vn = separate_vertices(F,V)
Cn = simplex2vertexdata(Fn,C,Vn)

# Visualization

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A multi-object mesh")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax2 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Grouping")
hp2 = poly!(ax2,GeometryBasics.Mesh(Vn,Fn), strokewidth=1,color=Cn, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(:viridis))

Legend(fig[1, 2],[hp1,hp2],["Mesh object","Grouped mesh"])
Colorbar(fig[2, 2],hp2, label = "Group labelling")

fig


