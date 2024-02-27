using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra

# Example geometry
testCase = 7
if testCase == 1
    s=1.0
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)

    F = Vector{TriangleFace{Int64}}(undef,1)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    # F,V=subtri(F,V,2)
elseif testCase==2
    s=1.0
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,5)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,    s, 0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s, 0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(   s,  0.0, 0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}( 2*s,    s, 0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}( 2*s,    -s, 0.0)

    F = Vector{TriangleFace{Int64}}(undef,2)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    F[2 ] = TriangleFace{Int64}(3,4,5)
    # F,V=subtri(F,V,2)
elseif testCase==3
    r = 1.0
    F,V = geosphere(3,r)    
elseif testCase==4
    r = 1.0
    F,V = quadsphere(3,r)    
elseif testCase==5 # Unmerged STL, each triangle is seperate group
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
elseif testCase==6 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase==7 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
end

indStart = [1]
# con_V2V = con_vertex_vertex_f(F,V)
# d_V2V = Dict{Vector{Int64},Float64}() 
# for i ∈ eachindex(V)
#     for ii ∈ con_V2V[i]
#         k = sort([i,ii])
#         if !haskey(d_V2V,k)
#             d_V2V[sort([i,ii])] = norm(V[i]-V[ii])
#         end 
#     end
# end

# con_V2V = con_vertex_vertex_f(F,V)
# d_V2V = [ [norm(V[i]-V[ii]) for ii ∈ con_V2V[i]] for i ∈ eachindex(V)]

function distmarch(F,V,indStart; dist_tol=1e-3, numPoints=1)
    dd = Dict{Vector{Int64},Float64}() 
    d = fill(Inf,length(V))
    d[indStart] .= 0.0
    ind = Set{Int64}();
    for i ∈ indStart
        push!(ind,i)
    end
    indSeed = Vector{Int64}()
    con_V2V = con_vertex_vertex_f(F,V) # vertex face connectivity 
    c = false
    ds = -1.0 # Set negative initially 
    while true          
        for i ∈ ind
            for ii ∈ con_V2V[i]
                k = sort([i,ii])
                if !haskey(dd,k)                    
                    dd[k] = norm(V[i]-V[ii])
                    push!(ind,ii)   
                end
                d[ii] = min(d[ii],dd[k]+d[i])  
            end
        end
        if ~any(isinf.(d)) # Start checking once all are no longer Inf
            if c # If we were here before
                if abs(sum(d)-ds)<dist_tol                    
                    if length(indSeed)==numPoints
                        break
                    else 
                        indNewSeed = findmax(d)[2]
                        push!(indSeed,indNewSeed)
                        d[indNewSeed] = 0.0                                           
                    end
                end
            end
            c = true # Flip to denote we've been here           
            ds = sum(d) # Now start computing the sum to check convergence
        end
    end
    return d, indSeed
end

function pointseed(F,V,numPoints; indStart = 1, dist_tol=1e-3)    
    ind = [indStart]
    for _ = 1:1:numPoints
        d = distmarch(F,V,ind; dist_tol=dist_tol)
        push!(ind,findmax(d)[2])
    end
    return ind, distmarch(F,V,ind; dist_tol=dist_tol)
end

numPoints = 10
# ind,d = pointseed(F,V,numPoints; indStart = 1)

d,ind = distmarch(F,V,indStart; numPoints = numPoints)

# c = cgrad(Reverse(:Spectral),numPoints,categorical = true)

M=GeometryBasics.Mesh(V,F)
fig = Figure(size=(800,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Distance marching")


# hp1 =poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=0.5,color=d, strokecolor=:black, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
hp1 =mesh!(ax1,GeometryBasics.Mesh(V,F), color=d, shading = FastShading, transparency=true,colormap=Reverse(:Spectral))

scatter!(ax1,V[ind],color=:black,markersize=15)
# scatter!(ax1,V,color=d,markersize=10)

# scatter!(ax1,V[indNext],color=:blue,markersize=30)

# normalplot(ax1,M,color=:red,linewidth=3)
# ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Grouping")
# for q = 1:1:maximum(C)
#     poly!(ax2,GeometryBasics.Mesh(V,F[C.==q]), strokewidth=2,color=c[q], strokecolor=:black, shading = FastShading, transparency=false)
# end
Colorbar(fig[1, 2], hp1)
fig


