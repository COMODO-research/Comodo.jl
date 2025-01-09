using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `boundaryfaces` to obtain the boundary faces of a 
volumetric mesh. 
=#

testCase = 1

if testCase == 1 
    sampleSize = 10
    pointSpacing = 2
    boxDim = sampleSize.*[1,1,1] # Dimensionsions for the box in each direction
    boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
    E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)    
elseif testCase == 2
    F1,V1 = geosphere(3,1.0)
    E,V,CE,Fb,Cb = tetgenmesh(F1,V1)   
    F = element2faces(E)
end

# Get boundary faces from either the total set of element faces or the elements
Fb1 = boundaryfaces(F)
Fb2 = boundaryfaces(E)

# Visualisation

fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Input faces")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=true, overdraw=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary faces")
hp1 = poly!(ax2,GeometryBasics.Mesh(V,Fb2), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=true, overdraw=false)

fig