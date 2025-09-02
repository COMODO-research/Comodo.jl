using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# Example curve
r = 1.0
nc = 16
Vc = circlepoints(r,nc;dir=:acw)

d = 3.0 # Extrusion distance (extent)
n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0)) # Extrusion direction
direction = :positive

F1,V1 = extrudecurve(Vc; extent=d, direction=:positive, n=n, close_loop=true, face_type=:quad)
F2,V2 = extrudecurve(Vc; extent=d, direction=:both, n=n, close_loop=true, face_type=:forwardslash)
F3,V3 = extrudecurve(Vc; extent=d, direction=:negative, n=n, close_loop=true, face_type=:tri_even)

n = normalizevector(Vec{3, Float64}(1.0,0.0,1.0))
F4,V4 = extrudecurve(Vc; extent=d, direction=:positive, n=n, close_loop=true, face_type=:quad)
F5,V5 = extrudecurve(Vc; extent=d, direction=:both, n=n, close_loop=true, face_type=:tri_even)
F6,V6 = extrudecurve(Vc; extent=d, direction=:negative, n=n, close_loop=true, face_type=:quad2tri)

## Visualization
GLMakie.closeall()

fig = Figure(size=(1200,1200))

ax1 = AxisGeom(fig[1, 1], limits=(-r,r,-r,r,-d,d), title = """Extruded direction=:positive, face_type=:quad """)
hp1 = lines!(ax1,Vc,color=:red,linewidth=4, transparency=true, depth_shift=-1.0f-3)
hp2 = meshplot!(ax1, F1, V1)

ax2 = AxisGeom(fig[1, 2], limits=(-r,r,-r,r,-d,d), title = """Extruded direction=:both, face_type=:forwardslash """)
hp1 = lines!(ax2,Vc,color=:red,linewidth=4, transparency=true, depth_shift=-1.0f-3)
hp2 = meshplot!(ax2, F2, V2)

ax3 = AxisGeom(fig[1, 3], limits=(-r,r,-r,r,-d,d), title = """Extruded direction=:negative, face_type=:tri_even """)
hp1 = lines!(ax3,Vc,color=:red,linewidth=4, transparency=true, depth_shift=-1.0f-3)
hp2 = meshplot!(ax3, F3, V3)

ax4 = AxisGeom(fig[2, 1], title = """Extruded direction=:positive, face_type=:quad """)
hp1 = lines!(ax4,Vc,color=:red,linewidth=4, transparency=true, depth_shift=-1.0f-3)
hp2 = meshplot!(ax4, F4, V4)

ax5 = AxisGeom(fig[2, 2], aspect = :data,title = """Extruded direction=:both, face_type=:tri_even """)
hp1 = lines!(ax5,Vc,color=:red,linewidth=4, transparency=true, depth_shift=-1.0f-3)
hp2 = meshplot!(ax5, F5, V5)

ax6 = AxisGeom(fig[2, 3], title = """Extruded direction=:negative, face_type=:quad2tri """)
hp1 = lines!(ax6,Vc,color=:red,linewidth=4, transparency=true, depth_shift=-1.0f-3)
hp2 = meshplot!(ax6, F6, V6)

fig