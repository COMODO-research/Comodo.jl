using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

testCase = 1
if testCase == 1
    fileName_mesh = joinpath(comododir(),"assets","obj","spot_quadrangulated.obj")
    fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
else
    fileName_mesh = joinpath(comododir(),"assets","obj","lego_figure.obj")
    fileName_texture = joinpath(comododir(),"assets","obj","lego_figure.png")    
end

M = load(fileName_mesh)
F = tofaces(GeometryBasics.faces(M))
V = topoints(GeometryBasics.coordinates(M))
# M = GeometryBasics.Mesh(V,F)

## Visualization
fig = Figure(size=(800,800))

# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Spot the cow")
ax1 = LScene(fig[1,1])
hp1 = poly!(ax1,M, color=load(fileName_texture), shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
# normalplot(ax1,M)
fig