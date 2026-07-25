using Comodo
using Comodo.GLMakie

r = 1.0 # Cylinder radius
h = 3.0 # Cylinder height

pointSpacing = 0.5 # Point spacing to use 
nr = ceil(Int, (2π*r)./pointSpacing) # Derived number of radial points to approximate point spacing
nh = 1 + ceil(Int, h./pointSpacing) # Derived number of height direction points to approximate point spacing

F1, V1 = cylinder(r, h, nr, nh; direction=:positive, face_type=:quad, face_orientation=:outward)
F2, V2 = cylinder(r, h, nr, nh; direction=:both, face_type=:forwardslash, face_orientation=:outward)
F3, V3 = cylinder(r, h, nr, nh; direction=:negative, face_type=:tri, face_orientation=:inward)

nh4 = 1 + ceil(Int, h./(pointSpacing*sqrt(3.0)/2.0))
F4, V4 = cylinder(r, h, nr, nh4; direction=:both, face_type=:backslash, face_orientation=:outward)

## Visualization
GLMakie.closeall()

fig = Figure(size=(1200,800))

ax1 = AxisGeom(fig[1, 1], title = """direction=:positive, face_type=:quad, face_orientation=:outward""")
hp1 = meshplot!(ax1, F1, V1)
normalplot(ax1, F1, V1)

ax2 = AxisGeom(fig[1, 2], title = """direction=:both, face_type=:forwardslash, face_orientation=:outward""")
hp1 = meshplot!(ax2, F2, V2)
normalplot(ax2, F2, V2)

ax3 = AxisGeom(fig[2, 1], title = """direction=:negative, face_type=:tri, face_orientation=:inward""")
hp1 = meshplot!(ax3, F3, V3)
normalplot(ax3, F3, V3)

ax4 = AxisGeom(fig[2, 2], title = """direction=:negative, face_type=:backslash, face_orientation=:outward""")
hp1 = meshplot!(ax4, F4, V4)
normalplot(ax4, F4, V4)

fig