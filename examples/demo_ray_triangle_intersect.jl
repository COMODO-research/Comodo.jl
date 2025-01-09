using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
# using FileIO

# Example geometry
F,V = geosphere(1,1.0)

M = GeometryBasics.Mesh(V,F)

## Visualization

np = 11
global markerSize = 20

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :ray, triSide=1""")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end

ax1 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :ray, triSide=0""")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 0)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end

ax1 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :ray, triSide=-1""")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = -1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end



ax1 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :line, triSide=1""")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end

ax1 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :line, triSide=0""")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 0)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end

ax1 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :line, triSide=-1""")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = -1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:red)
    scatter!(ax1,P,markersize = markerSize,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end


fig
