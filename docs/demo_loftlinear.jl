using Gibbon
using GLMakie
using GeometryBasics
using LinearAlgebra

## Create example data
r = 1.0 # Radius
nc = 35 # Number of points

V1 = circlepoints(r,nc)
height = 2.0
rFun(t) = r+0.5*sin(3*t)
V2 = [GeometryBasics.Point{3, Float64}(rFun(t)*cos(t),rFun(t)*sin(t),height) for t ∈ range(0,2*π-(2*π)/nc,nc)]
## Loft from curve 1 to curve 2
num_loft = 12
close_loop = true

F,V = loftlinear(V1,V2;num_loft=num_loft,close_loop=close_loop)

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A lofted surface")
hp1=lines!(ax1,V1, linewidth = 6, color = :blue)
scatter!(V1,markersize=25,color = :blue)
hp2=lines!(ax1,V2, linewidth = 6, color = :red)
scatter!(V2,markersize=25,color = :red)
hp3=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white,shading=FastShading,transparency=false)
# normalplot(ax1,GeometryBasics.Mesh(V,F); type_flag="face", color=:black)
Legend(fig[:, 2],[hp1,hp2,hp3],["curve 1", "curve 2", "lofted surface"])

fig