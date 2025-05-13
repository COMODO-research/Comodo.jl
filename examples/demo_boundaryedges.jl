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

# Using `boundaryedges` to find the boundary edges (edges only touching one face)
Eb = boundaryedges(F) # or equivalently Eb = boundaryedges(meshedges(F))

ind = edges2curve(Eb)
## Visualization
fig = Figure(size=(1200,800))
ax1 = AxisGeom(fig[1, 1]; title = "Boundary curve")
hp1 = meshplot!(ax1, F, V)
hp2 = edgeplot!(ax1, Eb, V; color=:red, linewidth=4.0f0)
Legend(fig[1, 2],[hp1,hp2],["Surface","Boundary edges"])
fig
