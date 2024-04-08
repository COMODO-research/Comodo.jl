using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

## Create example data
r = 1.0 # Radius
nc = 26 # Number of points

V1 = circlepoints(r,nc;dir=:cw)
height = 2.5
rFun(t) = r + 0.5.* sin(4.0*t)
V2 = circlepoints(rFun,nc;dir=:cw)
V2 = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V2]
V2,_ = evenly_sample(V2,nc)
## Loft from curve 1 to curve 2
num_steps = 11 # Uneven works best for "tri" face type
close_loop = true


face_types=[:quad,:tri_slash,:tri,:quad2tri]

## Visualization
fig = Figure(size=(1200,600))

for q = 1:4
    F,V = loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type=face_types[q])

    ax1 = Axis3(fig[1, q], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A lofted surface: " * string(face_types[q]))
    hp1=lines!(ax1,V1, linewidth = 6, color = :blue)
    scatter!(V1,markersize=25,color = :blue)
    hp2=lines!(ax1,V2, linewidth = 6, color = :red)
    scatter!(V2,markersize=25,color = :red)
    hp3=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white,shading=FastShading,transparency=false)
    # normalplot(ax1,GeometryBasics.Mesh(V,F); type_flag=:face, color=:black)
end
# Legend(fig[1, 4],[hp1,hp2,hp3],["curve 1", "curve 2", "lofted surface"])

fig