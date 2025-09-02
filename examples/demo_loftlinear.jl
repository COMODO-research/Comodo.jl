using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

## Create example data
r = 1.0 # Radius
nc = 26 # Number of points

V1 = circlepoints(r,nc;dir=:acw)
height = 2.5
rFun(t) = r + 0.5.* sin(4.0*t)
V2 = circlepoints(rFun,nc;dir=:acw)
V2 = [GeometryBasics.Point{3, Float64}(v[1],v[2],height) for v ∈ V2]
V2 = evenly_sample(V2,nc)
## Loft from curve 1 to curve 2
num_steps = 11 # Uneven works best for "tri" face type
close_loop = true

face_types = [:quad,:forwardslash,:backslash,:tri,:tri_even,:quad2tri]

## Visualization
GLMakie.closeall()

markersize = 10
linewidth = 2 

fig = Figure(size=(1200,600))

nRows = 2
for (q,face_type) in enumerate(face_types)
    i = mod1(q,nRows)
    j = ceil.(Int,q/nRows)

    F,V = loftlinear(V1,V2; num_steps=num_steps, close_loop=close_loop, face_type=face_type)

    ax1 = AxisGeom(fig[i, j], title = "Lofted surface, face_type=:$face_type")    
    hp1 = lines!(ax1,V1, linewidth = linewidth, color = :blue)
    hp2 = scatter!(V1,markersize=markersize,color = :blue)
    hp3 = lines!(ax1,V2, linewidth = linewidth, color = :red)
    hp4 = scatter!(V2,markersize=markersize,color = :red)
    hp5 = meshplot!(ax1,F,V)
    # normalplot(ax1,GeometryBasics.Mesh(V,F); type_flag=:face, color=:black)
end

fig