using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `dist` to compute distances for ND points. A 3D 
point set is defined on an icosahedron. Next a refined (using subtri) version is
created and the distance from the refined to the unrefined are computed. Next
the minimum distances are visualised on the mesh. 
=#

# Defining icosahedron
r = 1 # radius of icosahedron
n = 2 # Number of refinement steps

# Define an icosahedron
F,V = platonicsolid(4,r) # GeometryBasics mesh description of icosahedron

# Created refined version
Fn,Vn = subtri(F,V,n) # Subdevide/refine the mesh linearly 
F,V = subtri(F,V,1)
# Use distND to compute distances from all in set 1 to all in set 2
DD = dist(Vn,V)

# Computed minimal distances. Note this is equivalent to: Dn = minDist(Vn,V; getIndex = false)
Dn = minimum(DD,dims=2)[:,1] 

# Visualization
GLMakie.closeall()

fig = Figure(size = (800,800))

ax1 = AxisGeom(fig[1, 1])
hp = meshplot!(ax1, Fn, Vn, color=Dn, colormap=:Spectral)
hs1 = scatter!(ax1, V,markersize=35,color=:black)
hs2 = scatter!(ax1, Vn,markersize=15,color=Dn,colormap=:Spectral)

Colorbar(fig[1, 2], hp,label="Distance")
Legend(fig[1, 3],[hp,hs1],["Distances on mesh","Point set"])

ax2 = Axis(fig[2, 1], aspect = DataAspect(), xlabel = "Point indices set 1 (refined)",
        ylabel = "Point indices set 2(icosahedron)")
hi = image!(ax2, DD, colormap=:Spectral,interpolate=false)
Colorbar(fig[2, 2], hi,label="Distance")
fig