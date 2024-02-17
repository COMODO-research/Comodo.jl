using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using LinearAlgebra
using FileIO

# Example geometry
testCase = 1
if testCase == 1
    F,V = geosphere(5,1.0)
elseif testCase == 2
    # Loading a mesh
    fileName_mesh = joinpath(gibbondir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = faces(M)
    V = coordinates(M)
    F = togeometrybasics_faces(F)
    V = togeometrybasics_points(V)
    F,V,ind1,ind2 = mergevertices(F,V)
    # F,V=subtri(F,V,1)
elseif testCase == 3
    # Loading a mesh
    fileName_mesh = joinpath(gibbondir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = faces(M)
    V = coordinates(M)
    F = togeometrybasics_faces(F)
    V = togeometrybasics_points(V)
    F,V,ind1,ind2 = mergevertices(F,V)
end

M = GeometryBasics.Mesh(V,F)


function ray_triangle_intersect(f::TriangleFace{Int64},V,ray_origin,ray_vector; rayType = "ray", triSide = 1, tolEps = eps(Float64))
    """
    Implementation of the Möller-Trumbore triangle-ray intersection algorithm. 

    Möller, Tomas; Trumbore, Ben (1997). "Fast, Minimum Storage Ray-Triangle Intersection".
    Journal of Graphics Tools. 2: 21-28. doi:10.1080/10867651.1997.10487468.    
    """
    
    # Edge vectors
    P1 = V[f[1]] # First corner point
    vec_edge_1 = V[f[2]].-P1 # Edge vector 1-2
    vec_edge_2 = V[f[3]].-P1 # Edge vector 1-3

    # Determine if ray/lines is capable of intersecting based on direction
    ray_cross_e2 = cross(ray_vector,vec_edge_2) 
    det_vec = dot(vec_edge_1,ray_cross_e2)  # Determinant det([-n' P21' P31'])
    if triSide == 1 # Pointing at face normals
        boolDet = det_vec>tolEps
    elseif triSide == 0 # Both ways
        boolDet = abs(det_vec)>tolEps
    elseif triSide == -1 # Pointing allong face normals
        boolDet = det_vec<tolEps
    end

    p = GeometryBasics.Point{3, Float64}(NaN,NaN,NaN)
    if boolDet        
        s = ray_origin.-P1
        u = dot(s,ray_cross_e2)/det_vec    
        if u >= 0 && u <= 1 # On triangle according to u            
            s_cross_e1 = cross(s,vec_edge_1)
            v = dot(ray_vector,s_cross_e1)/det_vec
            if v >= 0 && (u+v) <= 1 # On triangle according to both u and v
                # Allong ray/line coordinates i.e. intersection is at ray_origin + t.*ray_vector 
                t = dot(vec_edge_2,s_cross_e1)/det_vec                      
                if rayType == "ray" || (rayType == "line" && t>=0 && t<=1.0)                                                   
                    p = ray_origin .+ t.*ray_vector # same as: push!(P, P1 .+ u.*P21 .+ v.*P31)            
                end
            end
        end    
    end    
    return p 
end

function ray_triangle_intersect(F::Vector{TriangleFace{Int64}},V,ray_origin,ray_vector; rayType = "ray", triSide = 1, tolEps = eps(Float64))
    P = Vector{GeometryBasics.Point{3, Float64}}()
    indIntersect = Vector{Int64}()
    for qf ∈ eachindex(F)
        p = ray_triangle_intersect(F[qf],V,ray_origin,ray_vector; rayType = rayType, triSide = triSide, tolEps = tolEps)        
        if ~any(isnan.(p))
            push!(P,p)
            push!(indIntersect,qf)
        end
    end
    return P,indIntersect
end


## Visualization

fig = Figure(size=(800,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")
hp1 = poly!(ax1,M,color=:white, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.5)
# hp2 = normalplot(ax1,M,color=:red)

for x = range(-2,2,50)
    ray_origin = GeometryBasics.Point3{Float64}(x,0.0,2)
    ray_vector = Vec3{Float64}(0.0,0.0,-2.75)
    P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector;rayType = "ray", triSide = 1)
    scatter!(ax1,ray_origin,markersize=25,color=:blue)
    scatter!(ax1,ray_origin.+ray_vector,markersize=25,color=:red)
    scatter!(ax1,P,markersize=36,color=:green)
    lines!(ax1,[ray_origin,ray_origin.+ray_vector],color=:blue)
    poly!(ax1,GeometryBasics.Mesh(V,F[indIntersect]),shading = FastShading, transparency=false, color=:green,strokecolor=:green, strokewidth=2)
end

fig
