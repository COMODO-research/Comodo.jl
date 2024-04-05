using Comodo
using GLMakie
using GeometryBasics
using FileIO

fileName_mesh = joinpath(comododir(),"assets","obj","spot_quadrangulated.obj")
# M = load(fileName_mesh; facetype = GeometryBasics.QuadFace{GeometryBasics.GLIndex})
M = load(fileName_mesh)
F = togeometrybasics_faces(faces(M))
V = togeometrybasics_points(coordinates(M))
# M = GeometryBasics.Mesh(V,F)

# F,V = mergevertices(F,V)

# M=cube(1)


## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Spot the cow")
hp1 = poly!(ax1,M, color=load(joinpath(comododir(),"assets","obj","spot_texture.png")), shading = FastShading, transparency=false,strokecolor=:black,strokewidth=2)
# hp1 = mesh!(ax1,M, color=:red, shading = FastShading, transparency=false)
# hp1 = poly!(ax1,M, color=:red, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=2)
# normalplot(ax1,M)

fig