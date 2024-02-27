using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra
using ProgressMeter

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
    F,V = geosphere(4,r)    
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

function distmarch(F,V,indStart; d=missing, dd=missing, dist_tol=1e-3,con_V2V=missing,l=missing)

    # Get vertex-vertex connectivity
    if ismissing(con_V2V)
        con_V2V = con_vertex_vertex_f(F,V) 
    end

    # Compute "Laplacian umbrella" distances
    if ismissing(dd)
        dd = Dict{Vector{Int64},Float64}()  
        for i ∈ eachindex(V)
            for ii ∈ con_V2V[i]
                k = sort([i,ii])
                if !haskey(dd,k)
                    dd[sort(k)] = norm(V[i]-V[ii])
                end 
            end
        end
    end

    # Get/allocate distance vector
    if ismissing(d)
        d = fill(Inf,length(V))
    end

    if ismissing(l)
        l = fill(0,length(V))
    end

    # Set start distances to zero 
    d[indStart] .= 0.0
    l[indStart] .= 1:1:length(indStart)
    
    c = false
    ds = -1.0 # Set negative initially 

    boolCheck = fill(true,length(V))
    while true                          
        for i ∈ eachindex(V) # For each point            
            for ii ∈ con_V2V[i] # Check umbrella neighbourhood
                minVal,minInd = findmin([d[ii],dd[sort([i,ii])]+d[i]])            
                if minInd==2
                    d[ii] = minVal                          
                    l[ii] = l[i]
                end
            end            
        end
        if ~any(isinf.(d)) # Start checking once all are no longer Inf
            if c # If we were here before
                if abs(sum(d)-ds)<dist_tol                                        
                    break                    
                end
            end
            c = true # Flip to denote we've been here           
            ds = sum(d) # Now start computing the sum to check convergence
        end
    end
    return d,dd,l
end

function distseedpoints(F,V,numPoints; ind=[1],dist_tol=1e-3)
    
    con_V2V = con_vertex_vertex_f(F,V) 
    d,dd,l = distmarch(F,V,ind; dist_tol=dist_tol,con_V2V=con_V2V)

    if numPoints>1
        @showprogress 1 "<distseedpoints>: Seeding points..." for q ∈ 2:1:numPoints            
            push!(ind,findmax(d)[2])
            d,dd,l = distmarch(F,V,ind; dist_tol=dist_tol, dd=dd,d=d,con_V2V=con_V2V,l=l)        
        end
    end
    return ind,d,l
end

numPoints = 100
ind,d,l = distseedpoints(F,V,numPoints; ind=[1],dist_tol=1e-3)   

Vp = V[ind]
Fp = Vector{TriangleFace{Int64}}()
for f in F
    ii = l[f]
    if length(unique(ii))==3
        push!(Fp,TriangleFace{Int64}(ii))
    end
end
Mp=GeometryBasics.Mesh(Vp,Fp)

# c = cgrad(Reverse(:Spectral),numPoints,categorical = true)

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Distances")
# hp1 =poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=0.5,color=d, strokecolor=:black, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
hp1 =mesh!(ax1,GeometryBasics.Mesh(V,F), color=d, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
scatter!(ax1,V[ind],color=:black,markersize=15)
Colorbar(fig[1, 2], hp1)

ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Point regions")
# hp1 =poly!(ax2,GeometryBasics.Mesh(V,F), strokewidth=0.5,color=d, strokecolor=:black, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
hp1 =mesh!(ax2,GeometryBasics.Mesh(V,F), color=l, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
scatter!(ax2,V[ind],color=:black,markersize=15)
Colorbar(fig[1, 4], hp1)

ax2 = Axis3(fig[1, 5], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Reconstructed surface")
hp1 =poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth=0.5,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap=Reverse(:Spectral))
# scatter!(ax2,V[ind],color=:black,markersize=15)
# Colorbar(fig[1, 4], hp1)

fig


