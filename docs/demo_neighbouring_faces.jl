using GLMakie
using Gibbon
using GeometryBasics
using FileIO
using Statistics

"""
dfasdf
"""

function indTouch(E::Vector{Vector{Int}},F)
    A = Vector{Vector{Int64}}(undef,length(E))    
    for q ∈ eachindex(E)        
        A[q] = indTouch(E[q],F)       
    end
    return A
end

function indTouch(e::Vector{Int},F)        
    a = Vector{Int64}()   
    for qf ∈ eachindex(F)
        for qe ∈ eachindex(e)
            if e[qe] ∈ F[qf]                                            
                push!(a,qf)
                break
            end
        end
    end
    return a
end

function indTouch(e::Int,F)       
    a = Vector{Int64}()   
    for qf ∈ eachindex(F)        
        if e ∈ F[qf]            
            push!(a,qf)
        end        
    end
    return a
end

function mergeVertices(F,V; roundVertices = true, numDigitsMerge=missing)

    m = length(V)
    if roundVertices
        if ismissing(numDigitsMerge)
            E = meshEdges(F)
            d = [sqrt( sum((V[e[1]] .- V[e[2]]).^2) ) for e ∈ E]
            pointSpacing = mean(d)
            m = round(Int64,log10(pointSpacing))
            numDigitsMerge=6-m
        end

        # Create rounded coordinates to help obtain unique set
        VR = [round.(v,digits = numDigitsMerge) for v ∈ V]

        # Get unique indices and reverse for rounded vertices
        _,ind1,ind2 = gunique(VR; return_index=true, return_inverse=true,sort_entries=false)
        V = V[ind1] # The unique node set
    else
        V,ind1,ind2 = gunique(V; return_index=true, return_inverse=true,sort_entries=false)
    end

    if length(V) != m # If the length has changed
        # Correct indices for faces
        for q ∈ eachindex(F)
            F[q] = ind2[F[q]]
        end
    end

    return F,V,ind1,ind2
end

function faceEdgeFriends(F,indicesCheck = missing)

    if ismissing(indicesCheck)
        indicesCheck = 1:1:length(F)
    end
    
    A = Vector{Vector{Int64}}(undef,length(indicesCheck))
    for qi ∈ eachindex(indicesCheck)
        qff = indicesCheck[qi] 
        E = meshEdges(F[qff])
        a = Vector{Int64}(undef,length(E))   
        index_a = 1 
        for qf ∈ eachindex(F)            
            if qf!=qff                
                for e ∈ E
                    if all([in(i,F[qf]) for i ∈ e])
                        a[index_a] = qf
                        index_a += 1
                    end
                end
            end
            A[qi]= a 
        end
    end
    return A
end

# Loading a mesh
fileName_mesh = joinpath(gibbonDir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)
F = faces(M)
V = coordinates(M)
F,V,ind1,ind2 = mergeVertices(F,V; roundVertices=true)
V = toGeometryBasicsPoints(V)
F,V = subTri(F,V,1; method="loop")

M = GeometryBasics.Mesh(V,F)

# M = platonicsolid(4,1)
# F = faces(M)
# V = coordinates(M)
# F,V = subTri(F,V,2; method="loop")
# M = GeometryBasics.Mesh(V,F)

# Obtain mesh faces and vertices




indNeighbour = faceEdgeFriends(F)

## Visualisation

fig = Figure(size = (800,800))
ax = Axis3(fig[1, 1], aspect = :data)

hp1 = poly!(ax, M, strokewidth=2,color=:white,transparency=false,shading = FastShading)

poly!(ax, GeometryBasics.Mesh(V,F[1]), strokewidth=5,color=:green,transparency=true,strokecolor=:green,shading = FastShading)
poly!(ax, GeometryBasics.Mesh(V,F[indNeighbour[1]]), strokewidth=2,color=:red,transparency=true,shading = FastShading)

poly!(ax, GeometryBasics.Mesh(V,F[10]), strokewidth=5,color=:green,transparency=true,strokecolor=:green,shading = FastShading)
poly!(ax, GeometryBasics.Mesh(V,F[indNeighbour[10]]), strokewidth=2,color=:red,transparency=true,shading = FastShading)

# Legend(fig[1, 2],[hp1,hp2],["Initial","Rotated"])

fig


    