using GeometryBasics # For point and mesh format
using GLMakie
using Gibbon

"""
This demo shows the use of distND to compute distances for ND points. A 3D 
point set is defined on an icosahedron. Next a refined (using subtri) version is
created and the distance from the refined to the unrefined are computed. Next
the minimum distances are visualised on the mesh. 
"""

# Defining icosahedron
r = 1 # radius of icosahedron
n = 3 # Number of refinement steps

# Define an icosahedron
M = platonicsolid(4,r) # GeometryBasics mesh description of icosahedron
V = coordinates(M) # Get the mesh coordinates
F = faces(M) # Get the mesh faces

# Created refined version
Fn,Vn = subtri(F,V,n) # Subdevide/refine the mesh linearly 

# Compute nearest point distances
Dn = mindist(Vn,V; getIndex = false)

# Visualization
fig = Figure(size = (800,800))

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")

hp = poly!(ax,GeometryBasics.Mesh(Vn,Fn),strokewidth=0,color=Dn, shading=FastShading, overdraw=false)
hs1 = scatter!(ax, V,markersize=35,color=:black)
hs2 = scatter!(ax, Vn,markersize=15,color=Dn,strokewidth=1)

Colorbar(fig[1, 2], hs2,label="Distance")
Legend(fig[1, 3],[hp,hs1],["Distances on mesh","Point set"])
fig