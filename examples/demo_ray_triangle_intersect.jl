using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
# using FileIO

# Example geometry
F,V = geosphere(1, 1.0)
C = [sin(0.5*2π .* v[1]) for v in V] # Create example point data to interpolate later
# C = [1.0 for v in V] # Create example point data to interpolate later

M = GeometryBasics.Mesh(V,F)

function raytrace_interp(F, C, indIntersect, B)    
    CP = zeros(Float64, length(indIntersect))
    for (indexPoint, indexFace) in enumerate(indIntersect)                        
        for (i, j) in enumerate(collect(F[indexFace])) # For each node and corresponding barycentric coordinate/weight
            CP[indexPoint] += C[j] .* B[indexPoint][i] # Add contribution weighted by 
        end        
    end    
    return CP 
end


## Visualization
GLMakie.closeall()
cmap = :viridis
np = 11
global markerSize = 15

fig = Figure(size=(1200,800))

ax1 = AxisGeom(fig[1, 1], title = """rayType = :ray, triSide=1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=(:white,0.5), transparency=true)
# hp2 = normalplot(ax1,M,color=:red)

ray_vector = Vec3{Float64}(0.0, 0.0, -1.0)
for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)    
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 1)
    
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)    
    scatter!(ax1,P,markersize = markerSize,color=:red)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:red, strokecolor=:red, strokewidth=2)
end

ax1 = AxisGeom(fig[1, 2], title = """rayType = :ray, triSide=0""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=(:white,0.5), transparency=true)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)    
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 0)    

    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)    
    scatter!(ax1,P,markersize = markerSize,color=:red)
    meshplot!(ax1, F[indIntersect], V, color=:red, strokecolor=:red, strokewidth=2)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)    
end

ax1 = AxisGeom(fig[1, 3], title = """rayType = :ray, triSide=-1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=(:white,0.5), transparency=true)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)    
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = -1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)    
    scatter!(ax1,P,markersize = markerSize,color=:red)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:red, strokecolor=:red, strokewidth=2)
end


ax1 = AxisGeom(fig[1:2, 4], title = """rayType = :ray, triSide=1, barycentric interpolation""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=C, colormap=cmap, colorrange=[-1.0, 1.0], transparency=true) 
# hp2 = normalplot(ax1,M,color=:red)
Colorbar(fig[1:2, 5], hp1)
ray_vector = Vec3{Float64}(0.0, 0.0, -1.0)
for x = range(-1.25,1.25,np)
    ray_origin = Point3{Float64}(x,0.25*sin(x*pi),1.25)    
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 1)    
    
    scatter!(ax1, ray_origin, markersize = markerSize, color=:blue)       
    lines!(ax1,[ray_origin,ray_origin.+ray_vector], color=:blue)    

    if !isempty(P)
        CP = raytrace_interp(F, C, indIntersect, B)        
        println(CP)   
        scatter!(ax1, P, markersize = markerSize, color=CP, colormap=cmap, colorrange=[-1.0, 1.0], strokewidth=0.0)    
        meshplot!(ax1, F[indIntersect], V, color=(:red, 0.1), strokecolor=:red, strokewidth=2, transparency=true)
    end
end

ax1 = AxisGeom(fig[2, 1], title = """rayType = :line, triSide=1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=(:white,0.5), transparency=true)
# hp2 = normalplot(ax1,M,color=:red)
ray_vector = Vec3{Float64}(0.0, 0.0, -2.0)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)    
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 1)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:blue)
    scatter!(ax1,P,markersize = markerSize,color=:red)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:red, strokecolor=:red, strokewidth=2)
end

ax1 = AxisGeom(fig[2, 2], title = """rayType = :line, triSide=0""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=(:white,0.5), transparency=true)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 0)
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:blue)
    scatter!(ax1,P,markersize = markerSize,color=:red)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:red, strokecolor=:red, strokewidth=2)
end

ax1 = AxisGeom(fig[2, 3], title = """rayType = :line, triSide=-1""")
hp1 = meshplot!(ax1, F, V, strokewidth=0.5, color=(:white,0.5), transparency=true)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-1.25,1.25,np)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.25*sin(x*pi),1.25)
    ray_vector = Vec3{Float64}(0.0,0.0,-2)
    P, indIntersect, T, det_vals, B = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = -1)    
    scatter!(ax1,ray_origin,markersize = markerSize,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize = markerSize,color=:blue)
    scatter!(ax1,P,markersize = markerSize,color=:red)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    meshplot!(ax1, F[indIntersect], V, color=:red, strokecolor=:red, strokewidth=2)
end

fig