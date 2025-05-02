using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Statistics
using Comodo.GLMakie.Colors

#=
This demo shows the use of simplex2vertexdata to sample data for elements or 
faces onto the nodes instead. The sample features averaging of connected 
simplices and hence has a smoothing effect.  
=#

# Defining example mesh
r = 1.0 # radius
pointSpacing = 0.1
E,V = hexsphere(r,pointSpacing) # Volumetric elements 
F = boundaryfaces(E) # Surface elements 
L = meshedges(F)

# Create some example simplex data 
function exampleSimplexData(F,V)
    X = [v[1] for v in V]
    XF = [mean(X[f]) for f in F]  
    return sin.(2*2.0*pi*XF./r)
end
E_data = exampleSimplexData(E,V)
F_data = exampleSimplexData(F,V)
L_data = exampleSimplexData(L,V)

# Now use `simplex2vertexdata` to compute equivalent vertex data 
weighting = :size
E_data_on_V = simplex2vertexdata(E,E_data,V; weighting=weighting)
F_data_on_V = simplex2vertexdata(F,F_data,V; weighting=weighting)
E_data_on_V = simplex2vertexdata(L,L_data,V; weighting=weighting)

# Visualization
lineWidth = 2 
markersize = 15 
cAlpha = RGBA(1.0, 1.0, 1.0,0.25)
FE = element2faces(E)
FE_data = repeat(E_data,inner=6) 

E_s,VE_s = separate_vertices(E,V; scaleFactor = 0.75)
FE_s = element2faces(E_s)
FE_data_s = simplex2vertexdata(FE_s,FE_data)

F_s,VF_s = separate_vertices(F,V)
F_data_s = simplex2vertexdata(F_s,F_data)

L_s,VL_s = separate_vertices(L,V)
L_data_s = simplex2vertexdata(L_s,L_data)


fig = Figure(size = (1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title="Edge data")
hp1 = wireframe!(ax1,GeometryBasics.Mesh(VL_s,L_s), linewidth=lineWidth, color=L_data_s, colormap=:avocado)
Colorbar(fig[1, 2], hp1,label="Edge data")

ax2 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title="Node data")
mesh!(ax2,GeometryBasics.Mesh(V,F),color=cAlpha,shading=FastShading, overdraw=false,transparency=true)
hp2 = scatter!(ax2,V, color=E_data_on_V, colormap=:avocado, markersize=markersize)
Colorbar(fig[2, 2], hp2,label="Edge data on nodes")

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title="Face data")
hp3 = poly!(ax3,GeometryBasics.Mesh(VF_s,F_s), strokewidth=lineWidth,color=F_data_s,shading=FastShading, overdraw=false,colormap=:avocado)
Colorbar(fig[1, 4], hp3,label="Face data")

ax4 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title="Node data")
mesh!(ax4,GeometryBasics.Mesh(V,F),color=cAlpha,shading=FastShading, overdraw=false,transparency=true)
hp4 = scatter!(ax4,V, color=F_data_on_V, colormap=:avocado, markersize=markersize)
Colorbar(fig[2, 4], hp4,label="Face data on nodes")

ax4 = Axis3(fig[1, 5], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title="Element data")
hp5 = poly!(ax4,GeometryBasics.Mesh(VE_s,FE_s), strokewidth=0,color=FE_data_s,shading=FastShading, overdraw=false,colormap=:avocado)
Colorbar(fig[1, 6], hp3,label="Face data")

ax5 = Axis3(fig[2, 5], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title="Node data")
mesh!(ax5,GeometryBasics.Mesh(V,F),color=cAlpha,shading=FastShading, overdraw=false,transparency=true)
hp6 = scatter!(ax5,V, color=E_data_on_V, colormap=:avocado, markersize=markersize)
Colorbar(fig[2, 6], hp4,label="Face data on nodes")

fig