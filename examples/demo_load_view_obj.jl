using Comodo
using GLMakie
using GeometryBasics
using FileIO

fileName_mesh = joinpath(comododir(),"assets","obj","spot_quadrangulated.obj")
fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")

# M = load(fileName_mesh; facetype = GeometryBasics.QuadFace{GeometryBasics.GLIndex})

M = load(fileName_mesh)
F = tofaces(faces(M))
V = topoints(coordinates(M))
# M = GeometryBasics.Mesh(V,F)

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Spot the cow")
hp1 = poly!(ax1,M, color=load(fileName_texture), shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
# hp1 = mesh!(ax1,M, color=:red, shading = FastShading, transparency=false)
# hp1 = poly!(ax1,M, color=:red, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=2)
# normalplot(ax1,M)
fig