using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Statistics

#=
This demo shows the use of the `dirplot` function to visualize directional data. 
=#

fig = Figure(size=(1600,800))

F,V = icosahedron()
NV = vertexnormal(F,V; weighting=:area)

styleSet = (:to,:from,:through)

for i in eachindex(styleSet)
    ax1 = Axis3(fig[1, i], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = string(styleSet[i]))
    hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,shading=FastShading,color=:white, transparency=true, overdraw=false)        
    hpa = dirplot(ax1,V,NV.*mean(edgelengths(F,V)); color=:blue,linewidth=3,scaleval=1.0,style=styleSet[i])
end

fig
