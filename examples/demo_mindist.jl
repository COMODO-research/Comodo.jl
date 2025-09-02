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
Dn,indMin = mindist(Vn,V; getIndex = Val(true))

# Visualization
GLMakie.closeall()

fig = Figure(size = (1200,500))

ax1 = AxisGeom(fig[1, 1])

hp1 = meshplot!(ax1, Fn, Vn,strokewidth=0.0, color=Dn)
hs1 = scatter!(ax1, V,markersize=35,color=:black)

Colorbar(fig[1, 2], hp1,label="Distance")
Legend(fig[1, 3],[hp1,hs1],["Distances on mesh","Point set"])

ax2 = AxisGeom(fig[1, 4])

hp2 = meshplot!(ax2, Fn, Vn, strokewidth=0.0, color=indMin)
hs2 = scatter!(ax2, V,markersize=35,color=:black)

Colorbar(fig[1, 5], hp2,label="Nearest point index")
Legend(fig[1, 6],[hp2,hs2],["Point indices","Point set"])

fig