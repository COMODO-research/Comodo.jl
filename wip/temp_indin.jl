
function indIn(E::Vector{Vector{Int}},F)
    A = Vector{Vector{Int64}}(undef,length(E))    
    for q ∈ eachindex(E)        
        A[q] = indIn(E[q],F)       
    end
    return A
end

function indIn(e::Vector{Int},F)    
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

function indIn(e::Int,F)       
    indIn = Vector{Int64}()   
    for qf ∈ eachindex(F)        
        if e ∈ F[qf]            
            push!(indIn,qf)
        end        
    end
    return indIn
end