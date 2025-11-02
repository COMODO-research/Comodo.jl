using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

ns = 6
V1 = circlepoints(1.0, ns)
V2 = circlepoints(0.5, ns)
F, V = loftlinear(V1, V2; num_steps=nothing, close_loop=true, face_type=:quad)
F = [QuadFace{Int}(circshift(f,-1)) for f in F]

n = 1
method =:linear
direction1 = 0
Fn1,Vn1 = subquad(F, V, n; method=method, constrain_boundary=false, direction=direction1)  

direction2 = 1
Fn2,Vn2 = subquad(F, V, n; method=method, constrain_boundary=false, direction=direction2) 

direction3 = 2
Fn3,Vn3 = subquad(F, V, n; method=method, constrain_boundary=false, direction=direction3) 

## Visualization    
lineWidth = 4
fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Linear, n=$n, direction=$direction1")
hp1 = edgeplot!(ax1, F, V, linewidth=lineWidth, color=:red)
hp2 = meshplot!(ax1, Fn1, Vn1)

ax2 = AxisGeom(fig[1, 2], title = "Linear, n=$n, direction=$direction2")
hp3 = edgeplot!(ax2, F, V, linewidth=lineWidth, color=:red)
hp4 = meshplot!(ax2, Fn2, Vn2)

ax3 = AxisGeom(fig[1, 3], title = "Linear, n=$n, direction=$direction3")
hp5 = edgeplot!(ax3, F, V, linewidth=lineWidth, color=:red)
hp6 = meshplot!(ax3, Fn3, Vn3)

Legend(fig[1, 4],[hp1,hp2],["Initial","Refined"])

screen = display(GLMakie.Screen(), fig)