using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Statistics

#=
This demo shows the use of simplex2vertexdata to sample data for elements or 
faces onto the nodes instead. The sample features averaging of connected 
simplices and hence has a smoothing effect.  
=#

# Defining example mesh
r = 1 # radius
n = 3 # Number of refinement steps
F,V = geosphere(n,r) 

# Create face data 
X = [v[1] for v in V]
XF = [mean(X[f]) for f in F]  
DF = sin.(2*2.0*pi*XF./r)

DV = simplex2vertexdata(F,DF)

# Visualization
Fs,Vs = separate_vertices(F,V)
Cs = simplex2vertexdata(Fs,DF)

fig = Figure(size = (1200,500))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vs,Fs),strokewidth=0,color=Cs,
shading=FastShading, overdraw=false,colormap=:avocado)
Colorbar(fig[1, 2], hp1,label="Face data")

ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
hp2 = poly!(ax2,GeometryBasics.Mesh(V,F),strokewidth=0,color=DV,
shading=FastShading, overdraw=false,colormap=:avocado)
Colorbar(fig[1, 4], hp2,label="Vertex data")

fig