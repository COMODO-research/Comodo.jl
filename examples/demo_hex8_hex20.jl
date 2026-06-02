using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

pointSpacing = 0.5
boxDim = [2.5,3.1,4] # Dimensionsions for the box in each direction
boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 

E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

E_hex20, V_hex20 = hex8_hex20(E,V)

F_hex20 = element2faces(E_hex20)
indBoundary = boundaryfaceindices(F_hex20)

# Visualisation
cmap = Makie.Categorical(:Spectral) 

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,CFb_type)
M = GeometryBasics.Mesh(Vbs,Fbs)

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Boundary faces with boundary markers for the hexahedral mesh")
hp1 = meshplot!(ax1, Fb, V)
scatter!(ax1, V, color=:black, markersize=15)

ax2 = AxisGeom(fig[1, 2], title = "Boundary faces with boundary markers for the hexahedral mesh")
hp2 = meshplot!(ax2, F_hex20[indBoundary], V_hex20)
scatter!(ax2, V_hex20, color=:black, markersize=15)

fig