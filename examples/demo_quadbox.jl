using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

pointSpacing = 0.5
boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
F,V,C = quadbox(boxDim,boxEl)

## Visualization
GLMakie.closeall()

strokewidth1 = 2

cmap = Makie.Categorical(:Spectral) 

Fs,Vs = separate_vertices(F,V)
CV = simplex2vertexdata(Fs,C)

fig = Figure(size=(1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "A quadrangulated box")
hp1 = meshplot!(ax1, Fs, Vs, strokewidth=strokewidth1, color=CV, colormap=cmap)
# normalplot(ax1,Fs,Vs)
Colorbar(fig[1, 1][1, 2], hp1)
fig