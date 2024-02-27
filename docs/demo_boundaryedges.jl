using Comodo
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


# function edges2curve(Eb)
#     con_E2E = con_edge_edge(Eb)
#     ind = convert(Vector{Int64},Eb[1])
#     seen = fill(false,length(Eb))
#     i = 1 # Start with first edge
#     while !all(seen) # loop until all edges have been visited
#         ij = con_E2E[i]
#         seen[i] = true
#         if Eb[ij[1]][1]==ind[end]
#             i = ij[1]
#             push!(ind,Eb[i][2])
#             seen[i] = true
#         elseif Eb[ij[2]][1]==ind[end]
#             i = ij[2]
#             push!(ind,Eb[i][2])
#             seen[i] = true
#         end
#     end
#     return ind
# end
    # ind[end]
ind = edges2curve(Eb)

M = GeometryBasics.Mesh(V,F)

## Visualization
fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary curve Visualization")

hp1 = poly!(ax1,M, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
# hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,Eb), linewidth=5,color=:red)
hp2 = lines!(V[ind],color=:red,linewidth=5)
Legend(fig[1, 2],[hp1,hp2],["Surface","Boundary edges"])

fig
