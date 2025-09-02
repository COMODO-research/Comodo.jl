using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `subhex` to refine a a hexhedral mesh through splitting. 
=#

boxDim = [2.0,3.0,2.5] # Dimensions for the box in each direction
boxEl = [2,3,1] # Number of elements to use in each direction 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

Eh0,Vh0 = subhex(E,V,1;direction=0)
Fh0 = element2faces(Eh0)

Eh1,Vh1 = subhex(E,V,1;direction=1)
Fh1 = element2faces(Eh1)

Eh2,Vh2 = subhex(E,V,1;direction=2)
Fh2 = element2faces(Eh2)

Eh3,Vh3 = subhex(E,V,1;direction=3)
Fh3 = element2faces(Eh3)

# Visualisation
GLMakie.closeall()

strokewidth = 3
linewidth = 4

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Refined hexahedral mesh, direction=0 (all)")
hp1 = meshplot!(ax1, Fh0, Vh0, strokewidth=strokewidth)
hp2 = edgeplot!(ax1, F, V, linewidth=linewidth, color=:red, depth_shift=-0.01f0)
hpa = normalplot(ax1, Fh0, Vh0; color=:blue,linewidth=3)

ax2 = AxisGeom(fig[1, 2], title = "Refined hexahedral mesh, direction=1")
hp1 = meshplot!(ax2, Fh1, Vh1, strokewidth=strokewidth)
hp2 = edgeplot!(ax2, F, V, linewidth=linewidth, color=:red, depth_shift=-0.01f0)
hpa = normalplot(ax2, Fh1, Vh1; color=:blue,linewidth=3)

ax3 = AxisGeom(fig[2, 1], title = "Refined hexahedral mesh, direction=2")
hp1 = meshplot!(ax3, Fh2, Vh2, strokewidth=strokewidth)
hp2 = edgeplot!(ax3, F, V, linewidth=linewidth, color=:red, depth_shift=-0.01f0)
hpa = normalplot(ax3, Fh2, Vh2; color=:blue,linewidth=3)

ax4 = AxisGeom(fig[2, 2], title = "Refined hexahedral mesh, direction=3")
hp1 = meshplot!(ax4, Fh3, Vh3, strokewidth=strokewidth)
hp2 = edgeplot!(ax4, F, V, linewidth=linewidth, color=:red, depth_shift=-0.01f0)
hpa = normalplot(ax4, Fh3, Vh3; color=:blue,linewidth=3)

fig