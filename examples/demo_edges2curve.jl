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

# Using `boundaryedges` to find the boundary edges (edges only touching one face)
Eb = boundaryedges(F)

# Using `edges2curve` to convert the set of boundary edges to a list of points defining a curve
ind = edges2curve(Eb)

## Visualization
fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary edges converted to curve")

hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
hp2 = lines!(V[ind],color=:red,linewidth=6)
Legend(fig[1, 2],[hp1,hp2],["Surface","Boundary curve"])

fig
