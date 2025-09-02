using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# Example geometry for a sphere that is cut so some edges are boundary edges
nSub = 3 # Number of refinement steps of the geodesic sphere
r = 1.0 # Sphere radius
F,V = geosphere(nSub,r) # Creating the faces and vertices of a full sphere
VC = simplexcenter(F,V) # Finding triangle centre coordinates
F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces using z of central coordinates
F,V = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed

# get mesh edges
E = meshedges(F)

# Get unique edges
Eu,indReverse = gunique(E; return_unique=Val(true), return_inverse=Val(true), sort_entries=true)

# Count number of faces touching edges
C = count_edge_face(F,Eu,indReverse)

## Visualization
GLMakie.closeall()

Eun,Vn = separate_vertices(Eu,V)
Cn = simplex2vertexdata(Eun,C,Vn)

fig = Figure(size=(1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "Counts for faces connected to edges")

hp1 = meshplot!(ax1, F, V)
hp2 = edgeplot!(ax1, Eun, Vn, linewidth=5, color=Cn, colormap= colormap=Makie.Categorical(:viridis))

Legend(fig[1, 2],[hp1,hp2],["Surface","Edges"])
Colorbar(fig[1, 3],hp2, label = "Number of faces touching edge")
fig