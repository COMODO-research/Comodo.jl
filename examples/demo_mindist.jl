using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of dist to compute distances for ND points. A 3D 
point set is defined on an geodesic sphere. Next another, more refined version
is created, and the distance from this sphere to the coarser one is computed. 
Next the minimum distances are visualised on the mesh. 
=#

# Defining icosahedron
r = 1 # radius of icosahedron
n = 5 # Number of refinement steps

# Define a mesh
F,V = geosphere(1,r) 

# Define a mesh
Fn,Vn = geosphere(n,r)  # Subdevide/refine the mesh linearly 

# Compute nearest point distances
Dn,indMin = mindist(Vn,V; getIndex = Val(true))

# Visualization
fig = Figure(size = (1200,500))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")

hp1 = poly!(ax1,GeometryBasics.Mesh(Vn,Fn),strokewidth=0,color=Dn,shading=FastShading, overdraw=false)
hs1 = scatter!(ax1, V,markersize=35,color=:black)

Colorbar(fig[1, 2], hp1,label="Distance")
Legend(fig[1, 3],[hp1,hs1],["Distances on mesh","Point set"])

ax2=Axis3(fig[1, 4], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")

hp2 = poly!(ax2,GeometryBasics.Mesh(Vn,Fn),strokewidth=0,color=indMin,shading=FastShading, overdraw=false)
hs2 = scatter!(ax2, V,markersize=35,color=:black)

Colorbar(fig[1, 5], hp2,label="Nearest point index")
Legend(fig[1, 6],[hp2,hs2],["Point indices","Point set"])

fig