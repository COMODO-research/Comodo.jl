using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra
using Statistics

# Example geometry
# Example geometry
testCase = 1
if testCase == 1
    r = 25.25
    F,V = geosphere(5,r)  
    # V = [GeometryBasics.Point{3, Float64}(v[1],v[2],3*v[3]) for v ∈ V]  
elseif testCase==2
    r = 1.0
    F,V = quadsphere(3,r)    
elseif testCase==3
    r=sqrt(3)
    M = cube(r)
    F=faces(M)
    V=coordinates(M)
    # F = quad2tri(F,V; convert_method = "angle")
elseif testCase==4 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(gibbondir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
    F,V = subtri(F,V,2; method = "loop")
elseif testCase==5 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(gibbondir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase==6
    # Loading a mesh
    fileName_mesh = joinpath(gibbondir(),"assets","obj","motherChild_5k.obj")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))

end


M = GeometryBasics.Mesh(V,F)

A = facearea(F,V)


# EDGE-VERTEX connectivity
E = meshedges(F)
E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)  

# EDGE-FACE connectivity
con_E2F = con_edge_face(F,E_uni)
con_F2E = con_face_edge(F,E_uni,indReverse)                 
con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)


VE = map(e-> V[e[1]] .- V[e[2]],E_uni)
# DE = map(ve-> sqrt(sum(ve.^2)),VE)
edgeVectorLength = norm.(VE)
NE = VE ./ edgeVectorLength # normalize vector
N = facenormal(F,V)
edgeVectorLength = edgeVectorLength./mean(edgeVectorLength)

B = [acos(dot(N[con_E2F[q][1]],N[con_E2F[q][2]])) for q ∈ eachindex(E_uni)]

T = Vector{Matrix{Float64}}()
for q ∈ eachindex(E_uni)
    t = Matrix{Float64}(undef,(3,3))
    for i = 1:3
        for j=1:i
            t[i,j] = NE[q][i]*NE[q][j]*B[q]*edgeVectorLength[q]
            t[j,i] = t[i,j]
        end
    end
    push!(T,t)    
end

C_min = Vector{Float64}(undef,length(F))
C_max = Vector{Float64}(undef,length(F))
U_max = Vector{Vec{3,Float64}}(undef,length(F))
U_min = Vector{Vec{3,Float64}}(undef,length(F))
for q ∈ eachindex(F)
    n = length(con_F2E[q])
    t = sum(T[con_F2E[q]])./n # Average of edge tensors
    d,u = eigen(t)

    indSort1 = sortperm(abs.(d)) # sort so zero is first followed by min, max
    indSort2 = sortperm(d[indSort1[2:3]])
    indKeep = indSort1[2:3]
    indKeep = indKeep[indSort2]

    C_min[q] = d[indKeep[1]]
    C_max[q] = d[indKeep[2]]
    U_max[q] = u[:,indKeep[1]]
    U_min[q] = u[:,indKeep[2]]
 
end

c = [v[3] for v ∈ V]

# cp = Vec3()

C_min_V = simplex2vertexdata(F,C_min,V)
C_max_V = simplex2vertexdata(F,C_max,V)

VF = simplexcenter(F,V)

## Visualization

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Distance marching")
hp1 =mesh!(ax1,M, color=C_max_V, shading = FastShading, transparency=false,colormap=:bluesreds,colorrange = maximum(abs.(C_max_V)).*0.25.*(-1,1))
# hp = arrows!(ax1,VF,U_max./15,color=:black,quality=6) 

ax1 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Distance marching")
hp2 = mesh!(ax1,M, color=C_min_V, shading = FastShading, transparency=false,colormap=:bluesreds,colorrange = maximum(abs.(C_min_V)).*0.25.*(-1,1))
# hp = arrows!(ax1,VF,U_min./15,color=:black,quality=6) 

fig


