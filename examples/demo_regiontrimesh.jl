using Comodo
using GLMakie
using DelaunayTriangulation

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

n = 100
V1 = batman(n)

V2 = deepcopy(V1)
Vc = circlepoints(0.25,25; dir=:cw)
V2 = append!(V2,Vc)
V2 = [Point{2,Float64}(v[1],v[2]) for v in V2]

TR_1 = triangulate(V1)
F1 =  [TriangleFace{Int64}(tr) for tr in TR_1.triangles]

constrained_segments = append!(collect(1:length(V1)),1)
TR_2 = triangulate(V1; boundary_nodes=constrained_segments)
F2 =  [TriangleFace{Int64}(tr) for tr in TR_2.triangles]

constrained_segments1 = append!(collect(1:length(V1)),1)
constrained_segments2 = append!(collect(1:length(Vc)),1).+length(V1)
constrained_segments = [[constrained_segments1], [constrained_segments2]]
TR_3 = triangulate(V2; boundary_nodes=constrained_segments)
F3 =  [TriangleFace{Int64}(tr) for tr in TR_3.triangles]
V3 = [Point{3,Float64}(v[1],v[2],0.0) for v in V2]

TR_4 = deepcopy(TR_3)
stats = refine!(TR_4; max_area=1e-3, min_angle=27.3)
F4 =  [TriangleFace{Int64}(tr) for tr in TR_4.triangles]
V4 =  TR_4.points
V4 = [Point{3,Float64}(v[1],v[2],0.0) for v in V4]

E = boundaryedges(F4)
indB = unique(reduce(vcat,E))
n = 25
λ = 0.5
V4 = smoothmesh_laplacian(F4,V4, n, λ; constrained_points=indB)

# Visualisation
fig = Figure(size=(1800,600))

ax1 = Axis3(fig[1, 1],aspect = :data,title="Unconstrained")
hp1 = lines!(ax1, V1,linewidth=4,color=:blue)
hp2 = scatter!(ax1, V1,markersize=8,color=:red)
hp3 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=0.5,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax2 = Axis3(fig[1, 2],aspect = :data,title="Constrained")
hp1 = lines!(ax2, V1,linewidth=4,color=:blue)
hp2 = scatter!(ax2, V1,markersize=8,color=:red)
hp3 = poly!(ax2,GeometryBasics.Mesh(V1,F2), strokewidth=0.5,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax3 = Axis3(fig[2, 1],aspect = :data,title="Constrained with hole")
hp11 = lines!(ax3, V3[constrained_segments1],linewidth=4,color=:blue)
hp12 = lines!(ax3, V3[constrained_segments2],linewidth=4,color=:blue)
hp2 = scatter!(ax3, V3,markersize=8,color=:red)
hp3 = poly!(ax3,GeometryBasics.Mesh(V3,F3), strokewidth=0.5,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax4 = Axis3(fig[2, 2],aspect = :data,title="Constrained with hole")
hp11 = lines!(ax4, V4[constrained_segments1],linewidth=4,color=:blue)
hp12 = lines!(ax4, V4[constrained_segments2],linewidth=4,color=:blue)
hp2 = scatter!(ax4, V4,markersize=8,color=:red)
hp3 = poly!(ax4,GeometryBasics.Mesh(V4,F4), strokewidth=0.5,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

fig
