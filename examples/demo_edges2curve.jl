using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics


testCase = 3

if testCase == 1 # Full sphere, no boundary
    # Example geometry for a sphere
    nSub = 2 # Number of refinement steps of the geodesic sphere
    r = 1.0 # Sphere radius
    F,V = geosphere(nSub,r)

    # Using `boundaryedges` to find the boundary edges (edges only touching one face)
    Eb = boundaryedges(F)
elseif testCase == 2 # Full closed curve ()
    # Example geometry for a hemisphere 
    nSub = 2 # Number of refinement steps of the geodesic sphere
    r = 1.0 # Sphere radius
    F,V = hemisphere(nSub,r)

    # Using `boundaryedges` to find the boundary edges (edges only touching one face)
    Eb = boundaryedges(F)
elseif testCase == 3 # Partial non-closed curve (cut hemisphere with bottom edges excluded)
    tol_level = 1e-3

    # Example geometry for a sphere that is cut so some edges are boundary edges
    nSub = 2 # Number of refinement steps of the geodesic sphere
    r = 1.0 # Sphere radius
    F,V = hemisphere(nSub,r) # Creating the faces and vertices of a full sphere
    VC = simplexcenter(F,V) # Finding triangle centre coordinates
    F = [F[i] for i in findall(map(v-> v[1]>0,VC))] # Remove some faces using z of central coordinates
    F,V = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed
    # invert_faces!(F)

    # Using `boundaryedges` to find the boundary edges (edges only touching one face)
    Eb = boundaryedges(F)
    
    VC = simplexcenter(Eb,V) # Finding triangle centre coordinates
    Eb = [Eb[i] for i in findall(map(v-> v[3]>tol_level,VC))] # Remove some faces using z of central coordinates
    
end 

# Using `edges2curve` to convert the set of boundary edges to a list of points defining a curve
ind = edges2curve(Eb)

## Visualization
linewidth = 10
markersize = 35

fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary edges converted to curve")

hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth = 1, color = :white, strokecolor = :black, shading = FastShading, transparency = false)

if !isempty(ind)
    hp2 = lines!(ax1,V[ind], color = :red, linewidth = linewidth)
    hp3 = scatter!(ax1,V[ind[1]], color = :blue, markersize = markersize)
    hp4 = scatter!(ax1,V[ind[end]], color = :green, markersize = markersize)

    # hp1 = wireframe!(ax1,GeometryBasics.Mesh(V,Eb), linewidth=linewidth, color = :blue)
    Legend(fig[1, 2], [hp1, hp2, hp3, hp4], ["Surface","Boundary curve", "Curve start", "Curve end"])
else
    Legend(fig[1, 2], [hp1], ["Surface"])
end
fig