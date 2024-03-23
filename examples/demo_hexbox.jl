using Comodo
using GLMakie

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

sampleSize = 10
pointSpacing = 2
boxDim = sampleSize.*[1,1,1] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction 
# boxDim = (1.0,2.0,π) # Box dimensions
# boxEl = (2,4,6) # Number of hexahedral elements in each direction 

E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

M_F = GeometryBasics.Mesh(V,F)
M_Fb = GeometryBasics.Mesh(V,Fb)


# Visualisation

fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Wireframe of hexahedral mesh")
hp1 = wireframe!(ax1,M_F, linewidth=3,color=:red, overdraw=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary faces of hexahedral mesh")
hp2 = poly!(ax2,M_Fb, strokewidth=3,shading=FastShading,strokecolor=:red, color=:white, transparency=true, overdraw=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary faces and normals of hexahedral mesh")
hp3 = poly!(ax3,M_Fb, strokewidth=3,shading=FastShading,strokecolor=:red, color=:white, transparency=false, overdraw=false)
hp3 = normalplot(ax3,M_Fb; type_flag=:face, color=:black,linewidth=3)

fig
