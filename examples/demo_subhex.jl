using Comodo
using GLMakie
using GeometryBasics

#=
This demo shows the use of `subhex` to refine a a hexhedral mesh through splitting. 
=#

boxDim = [2.0,3.0,1.0] # Dimensions for the box in each direction
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
strokewidth = 3
linewidth = 4

fig = Figure(size=(1600,800))

ax0 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined hexahedral mesh, direction=0 (all)")
hp1 = poly!(ax0,GeometryBasics.Mesh(Vh0,Fh0), strokewidth=strokewidth,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
hp2 = wireframe!(ax0,GeometryBasics.Mesh(V,F), linewidth=linewidth,color=:red, overdraw=false)
hp3 = normalplot(ax0,Fh0,Vh0;color=:black)

ax1 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined hexahedral mesh, direction=1")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vh1,Fh1), strokewidth=strokewidth,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,F), linewidth=linewidth,color=:red, overdraw=false)
hp3 = normalplot(ax1,Fh1,Vh1;color=:black)

ax2 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined hexahedral mesh, direction=2")
hp1 = poly!(ax2,GeometryBasics.Mesh(Vh2,Fh2), strokewidth=strokewidth,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
hp2 = wireframe!(ax2,GeometryBasics.Mesh(V,F), linewidth=linewidth,color=:red, overdraw=false)
hp3 = normalplot(ax2,Fh2,Vh2;color=:black)

ax3 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined hexahedral mesh, direction=3")
hp1 = poly!(ax3,GeometryBasics.Mesh(Vh3,Fh3), strokewidth=strokewidth,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
hp2 = wireframe!(ax3,GeometryBasics.Mesh(V,F), linewidth=linewidth,color=:red, overdraw=false)
hp3 = normalplot(ax3,Fh3,Vh3;color=:black)

fig