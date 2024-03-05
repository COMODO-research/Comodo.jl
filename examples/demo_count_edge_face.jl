using Comodo
using GLMakie
using GeometryBasics

# Example geometry for a sphere that is cut so some edges are boundary edges
nSub = 3 # Number of refinement steps of the geodesic sphere
r = 1.0 # Sphere radius
F,V = geosphere(nSub,r) # Creating the faces and vertices of a full sphere
VC = simplexcenter(F,V) # Finding triangle centre coordiantes
F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces using z of central coordinates
F,V = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed

# get mesh edges
E = meshedges(F)

# Get unique edges
Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

# Count number of faces touching edges
C = count_edge_face(F,Eu,indReverse)

## Visualization

Eun,Vn = seperate_vertices(Eu,V)
Cn = simplex2vertexdata(Eun,C,Vn)

fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary curve Visualization")

hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
hp2 = wireframe!(ax1,GeometryBasics.Mesh(Vn,Eun), linewidth=5,color=Cn,colormap= colormap=Makie.Categorical(:viridis))

Legend(fig[1, 2],[hp1,hp2],["Surface","Edges"])
Colorbar(fig[1, 3],hp2, label = "Number of faces touching edge")
fig
