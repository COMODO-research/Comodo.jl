using Gibbon
using TetGen
using GeometryBasics
using GeometryBasics: Mesh, QuadFace
using GLMakie

# Construct a cube out of Quads
points = Point{3, Float64}[
    (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
    (2.0, 2.0, 0.0), (0.0, 2.0, 0.0),
    (0.0, 0.0, 2.0), (2.0, 0.0, 2.0),
    (2.0, 2.0, 2.0), (0.0, 2.0, 2.0)
]

facets = QuadFace{Cint}[
    1:4,
    5:8,
    [1,5,6,2],
    [2,6,7,3],
    [3, 7, 8, 4],
    [4, 8, 5, 1]
]

markers = Cint[-1, -2, 0, 0, 0, 0]
# attach some additional information to our faces!
mesh = Mesh(points, meta(facets, markers=markers))
result = tetrahedralize(mesh, "vpq1.414a0.1")


#Visualize mesh
GLMakie.activate!(inline=false) # To avoid plotting in plotpane as per: https://github.com/MakieOrg/Makie.jl/issues/2956
fig = Figure()
ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tetrahedron")

mesh!(ax1,normal_mesh(result), color=:blue, transparency=false)
wireframe!(ax1,result)

fig