using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
pointSpacing = 0.5

F,V,C = tribox(boxDim,pointSpacing)

# Visualization
GLMakie.closeall()

strokewidth1 = 2
lineWidth = 4
cmap = Makie.Categorical(:Spectral) 

Fs,Vs = separate_vertices(F,V)
CV = simplex2vertexdata(Fs,C)

fig = Figure(size=(1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "A triangulated box")
hp1 = meshplot!(ax1, Fs, Vs, color=CV, colormap=cmap)
# normalplot(ax1,Fs,Vs)
Colorbar(fig[1, 2], hp1)
fig