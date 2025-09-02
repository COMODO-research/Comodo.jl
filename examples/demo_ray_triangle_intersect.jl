using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
# using FileIO

# Example geometry
F,V = geosphere(1,1.0)

M = GeometryBasics.Mesh(V,F)

## Visualization
GLMakie.closeall()

np = 11
global markerSize = 20

fig = Figure(size=(1200,800))

ax1 = AxisGeom(fig[1, 1], title = """rayType = :ray, triSide=1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:green, strokecolor=:green, strokewidth=2)
end

ax1 = AxisGeom(fig[1, 2], title = """rayType = :ray, triSide=0""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 0)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:green, strokecolor=:green, strokewidth=2)
end

ax1 = AxisGeom(fig[1, 3], title = """rayType = :ray, triSide=-1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = -1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:green, strokecolor=:green, strokewidth=2)
end



ax1 = AxisGeom(fig[2, 1], title = """rayType = :line, triSide=1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:green, strokecolor=:green, strokewidth=2)
end

ax1 = AxisGeom(fig[2, 2], title = """rayType = :line, triSide=0""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 0)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:green, strokecolor=:green, strokewidth=2)
end

ax1 = AxisGeom(fig[2, 3], title = """rayType = :line, triSide=-1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = -1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:green, strokecolor=:green, strokewidth=2)
end

fig