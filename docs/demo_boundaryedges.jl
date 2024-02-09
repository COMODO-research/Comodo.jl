using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra
using FileIO
using SparseArrays

# Example geometry
F,V = geosphere(3,1.0)
VC = simplexcenter(F,V)
F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces

F,V = remove_unused_vertices(F,V)
    
Eb = boundaryedges(F)

M = GeometryBasics.Mesh(V,F)

## Visualization
fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")

hp1 = poly!(ax1,M, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)
Legend(fig[1, 2],[hp1,hp2],["Surface","Boundary edges"])

fig
