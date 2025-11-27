using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO 

r = 10.0 # radius
n = 1 # Number of refinement iterations
F1,V1 = geosphere(n,r)
V1 = [v + Point{3,Float64}(2.0, 2.0, 2.0,) for v in V1]

# Loading a mesh
fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)
F2 = [TriangleFace{Int}(f) for f in faces(M)]
V2 = coordinates(M)
F2,V2 = mergevertices(F2,V2)

V1C = surface_centroid(F1, V1)
V2C = surface_centroid(F2, V2)

## Visualization parameters
marker_size = 15 

## Visualize mesh
GLMakie.closeall()

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Triangle centre points")
hp1 = meshplot!(ax1, F1, V1, color=(:white,0.5), transparency=true)
hs2 = scatter!(ax1, V1C, markersize=marker_size, color=:blue)

ax2 = AxisGeom(fig[1, 2], title = "Quad centre points")
hp2 = meshplot!(ax2, F2, V2, color=(:white,0.5), transparency=true)
hs2 = scatter!(ax2, V2C, markersize=marker_size, color=:blue)

fig