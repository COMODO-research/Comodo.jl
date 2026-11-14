using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.GLMakie.Colors


r1 = 3.0
r2 = 5.0
n1 = 16
n2 = 51 # Use an uneven number for :tri and :tri_even if a homogeneous result is needed
face_types = (:quad, :forwardslash, :backslash, :tri, :tri_even)
pointSpacings = (1.0,    0.25,   1.0,  0.5,   1.0) 
spacing_types = (:inner, :middle, :middle, :middle, :outer)

## Visualization
GLMakie.closeall()

fig1 = Figure(size=(1200,1200))

for (i, face_type) in enumerate(face_types)
    F, V = torus(r1, r2, n1, n2; face_type=face_type)

    ax = AxisGeom(fig1[mod1(i,2), mod1(i,3)], title = "face_type = :$face_type")
    meshplot!(ax, F, V, color = :white, strokewidth=1.0)
end

screen = display(GLMakie.Screen(), fig1)
GLMakie.set_title!(screen, "Number of points")

##

fig2 = Figure(size=(1200,1200))

for (i, pointSpacing) in enumerate(pointSpacings)

    spacing_type = spacing_types[i]
    F, V = torus(r1, r2, pointSpacing; face_type=:quad, spacing_type=spacing_type)

    ax = AxisGeom(fig2[mod1(i,2), mod1(i,3)], title = "pointSpacing = :$pointSpacing, spacing_type=:$spacing_type")
    meshplot!(ax, F, V, color = :white, strokewidth=1.0)
end

screen = display(GLMakie.Screen(), fig2)
GLMakie.set_title!(screen, "Point spacing")
