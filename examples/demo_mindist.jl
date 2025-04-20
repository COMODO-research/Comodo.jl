using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of mindist to compute nearest distances from one 
pointset to another. 
=#

# Defining icosahedron
r = 1 # radius of icosahedron
n = 5 # Number of refinement steps

# Define a mesh
F,V = geosphere(1,r) 

# Define a mesh
Fn,Vn = geosphere(n,r)  # Subdevide/refine the mesh linearly 

# Compute nearest point distances
Dn,indMin = mindist(Vn,V; getIndex = true)

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