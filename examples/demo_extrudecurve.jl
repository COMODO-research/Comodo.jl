using Comodo
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra

# Example curve
r = 1
nc = 16
Vc = circlepoints(r,nc;dir=:cw)

d = 3.0
n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))
pointSpacing =0.25 
s = 1

#   num_loft = ceil(Int64,d/pointSpacing)
F1,V1 = extrudecurve(Vc,d;s=1, n=n, close_loop=true,face_type=:quad)
F2,V2 = extrudecurve(Vc,d;s=0, n=n, close_loop=true,face_type=:tri)
F3,V3 = extrudecurve(Vc,d;s=-1, n=n, close_loop=true,face_type=:tri_slash)

M1 = GeometryBasics.Mesh(V1,F1)
M2 = GeometryBasics.Mesh(V2,F2)
M3 = GeometryBasics.Mesh(V3,F3)

## Visualization
fig = Figure(size=(1200,1200))

ax1 = Axis3(fig[1, 1], aspect = :data, limits=(-r,r,-r,r,-d,d),xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Extruded s=1, face_type="quad" """)
hp1 = lines!(ax1,Vc,color=:red,linewidth=6, transparency=true, depth_shift=-1.0f-3)
hp2 = poly!(ax1,M1, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, limits=(-r,r,-r,r,-d,d), xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Extruded s=0, face_type="tri" """)
hp1 = lines!(ax2,Vc,color=:red,linewidth=6, transparency=true, depth_shift=-1.0f-3)
hp2 = poly!(ax2,M2, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, limits=(-r,r,-r,r,-d,d), xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Extruded s=-1, face_type="tri_slash" """)
hp1 = lines!(ax3,Vc,color=:red,linewidth=6, transparency=true, depth_shift=-1.0f-3)
hp2 = poly!(ax3,M3, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

# Legend(fig[1, 4],[hp1,hp2],["Input curve","Surface"])

fig
