using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra
using ProgressMeter

# Example geometry
testCase = 6
if testCase == 1
    s=1.0
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)

    F = Vector{TriangleFace{Int}}(undef,1)
    F[1 ] = TriangleFace{Int}(1,2,3)
    # F,V=subtri(F,V,2)
elseif testCase==2
    # Bowtie mesh
    s=1.0
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}( 2*s,    s, 0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}( 2*s,    -s, 0.0)

    F = Vector{TriangleFace{Int}}(undef,2)
    F[1 ] = TriangleFace{Int}(1,2,3)
    F[2 ] = TriangleFace{Int}(5,4,3)
    # F,V=subtri(F,V,2)
elseif testCase==3
    r = 1.0
    F,V = geosphere(4,r)    
elseif testCase==4
    r = 1.0
    F,V = quadsphere(3,r)    
elseif testCase==5 # Unmerged STL, each triangle is seperate group
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
elseif testCase==6 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase==7 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_ = mergevertices(F,V)
end

z = [v[3] for v ∈ V]

ind=[findmin(z)[2]]

d,dd,l = distmarch(F,V,ind; dist_tol=1e-3)

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Distances")
hp1 = mesh!(ax1,GeometryBasics.Mesh(V,F), color=d, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
# scatter!(ax1,V[ind],color=:black,markersize=25)
Colorbar(fig[1, 2], hp1)

ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Point regions")
hp1 = mesh!(ax2,GeometryBasics.Mesh(V,F), color=l, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
scatter!(ax2,V[ind],color=:black,markersize=25)
Colorbar(fig[1, 4], hp1)

fig


