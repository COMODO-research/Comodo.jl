using Gibbon
using GLMakie
using GeometryBasics
using Statistics

## Define example input
r = 1.0 #radius
n = 3
F,V =geoSphere(n,r)
E = meshEdges(F)

 
function indIn2(E::Vector{Vector{Int}},F)
    A = Vector{Vector{Int64}}(undef,length(E))
    for q ∈ eachindex(E)
        A[q] = indIn(E[q],F)
    end
    return A
end

function indIn2(e,F)
    n = length(F[1])
    indIn = Vector{Int64}()
    for q ∈ eachindex(F)
        if length(union(e,F[q]))==n
            push!(indIn,q)
        end
    end
    return indIn
end

function indIn3(E::Vector{Vector{Int}},F)
    A = Vector{Vector{Int64}}(undef,length(E))
    for q ∈ eachindex(E)             
        A[q] = indIn3(E[q],F)        
    end
    return A
end

function indIn3(e,F)    
    indIn = Vector{Int64}()
    for q ∈ eachindex(F)
        if all([i ∈ F[q] for i in e])
            push!(indIn,q)
        end
    end
    return indIn
end

function indIn4(E::Vector{Vector{Int}},F)
    A = Vector{Vector{Int64}}(undef,length(E))    
    for q ∈ eachindex(E)        
        A[q] = indIn4(E[q],F)       
    end
    return A
end

function indIn4(e,F)    
    n = length(e)
    indIn = Vector{Int64}()   
    for qf ∈ eachindex(F)
        for qe ∈ eachindex(e)
            if e[qe] ∉ F[qf]
                break
            elseif qe==n
                push!(indIn,qf)
            end
        end
    end
    return indIn
end
