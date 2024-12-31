using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra
using Rotations

pointSpacing = 0.5
boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
F,V,C = quadbox(boxDim,boxEl)

## Visualization
strokewidth1 = 2
lineWidth = 4
cmap = Makie.Categorical(:Spectral) 

Fs,Vs = separate_vertices(F,V)
CV = simplex2vertexdata(Fs,C)

fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1][1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A quadrangulated box")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vs,Fs), strokewidth=strokewidth1,color=CV,shading=FastShading,transparency=false,colormap=cmap)
# normalplot(ax1,Fs,Vs)
Colorbar(fig[1, 1][1, 2], hp1)
fig