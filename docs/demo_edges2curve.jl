using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra
using FileIO
using SparseArrays

# Example geometry
F,V = geosphere(2,1.0)
VC = simplexcenter(F,V)
F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces
F,V = remove_unused_vertices(F,V)
    
Eb = boundaryedges(F) # Boundary edges

# Convert the set of boundary edges to a list of points defining a curve
ind = edges2curve(Eb)

M = GeometryBasics.Mesh(V,F)

## Visualization
fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary edge Visualization")

hp1 = poly!(ax1,M, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
# hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)
hp2 = lines!(V[ind],color=:red,linewidth=5, transparency=true, depth_shift=-1.0f-3)
Legend(fig[1, 2],[hp1,hp2],["Surface","Boundary edges"])

fig
