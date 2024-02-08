using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra
using FileIO
using SparseArrays

# Example geometry
testCase = 3
if testCase == 1
    F,V = geosphere(2,1.0)
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

# function trisurfslice(F,V; p=(0,0,0))
# Input parameters
p = mean(V,dims=1)[1]; # Point on cutting plane
n = normalizevector(Vec{3, Float64}(0.0,1.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

function trisurfslice1(F,V,n = (0.0,0.0,1.0), p = mean(V,dims=1); snapTolerance = 0, output_type="full")

    # Compute dot product with slicing vector/plane
    d = map(v-> dot(n,v.-p),V)
    if snapTolerance != 0.0
        d[abs.(d).<snapTolerance] .= 0
    end
    LV = d.<0 # Boolean for points under the plane

    if !all(LV) && any(LV)
        if output_type == "above"
            indKeep_F = findall([~any(LV[f]) for f ∈ F])
        elseif output_type == "below"
            indKeep_F = findall([all(LV[f]) for f ∈ F])
        elseif output_type == "full"
            indKeep_F = findall([all(LV[f]) || ~any(LV[f]) for f ∈ F])
        end

        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        con_E2F = con_edge_face(F,E_uni) # EDGE-FACE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse) # FACE-EDGE connectivity          

        numUnder_E = map(e-> sum(LV[e]),E_uni)         
        numOver_E = map(e-> sum(.~LV[e]),E_uni) 
        LE = numUnder_E.==1
        # println(typeof(LE))
        Ec = E_uni[LE]
        indF = reduce(vcat,con_E2F[LE])

        # Compute edge-plane intersection points
        Vp = Vector{Point{3,Float64}}(undef,length(Ec))
        for q ∈ eachindex(Ec)
            pe = V[Ec[q][1]] # Vector origin
            ne = V[Ec[q][2]].-pe # Edge vector
            de = dot(n,ne) # slope in n direction
            f = d[Ec[q][1]]/de
            Vp[q] = pe .- f.*ne
        end

        # Formulate edge to intersection point index vector
        indMap = zeros(Int64,length(E_uni))
        indMap[LE] = 1:1:length(Ec)
        indMap[LE] .+= length(V)

        # Append new points to V    
        Vn = append!(deepcopy(V),Vp) 
        
        Fn = F[indKeep_F]    
        for q ∈ indF
            f = F[q]
            ind_e = con_F2E[q]
            numUnder_e = numUnder_E[ind_e]
            ind2 = findfirst(numUnder_e.==2)
            if ~isnothing(ind2) # Found two under
                indNext1 = wrapindex(ind2+1,3)
                indNext2 = wrapindex(ind2+2,3)
                if output_type == "below" || output_type == "full"
                    fq = QuadFace{Int64}(f[ind2],f[indNext1],indMap[ind_e[indNext1]],indMap[ind_e[indNext2]])
                    ft1 = quad2tri([fq],Vn; convert_method = "angle")            
                    push!(Fn,ft1[1])
                    push!(Fn,ft1[2])
                end
                if output_type == "above" || output_type == "full"
                    ft2 = TriangleFace{Int64}(f[indNext2],indMap[ind_e[indNext2]],indMap[ind_e[indNext1]])
                    push!(Fn,ft2)
                end                        
            end
            numOver_e = numOver_E[ind_e]
            ind2 = findfirst(numOver_e.==2)
            if ~isnothing(ind2) # Found two over
                indNext1 = wrapindex(ind2+1,3)
                indNext2 = wrapindex(ind2+2,3)
                if output_type == "above" || output_type == "full"
                    fq = QuadFace{Int64}(f[ind2],f[indNext1],indMap[ind_e[indNext1]],indMap[ind_e[indNext2]])
                    ft1 = quad2tri([fq],Vn; convert_method = "angle")  
                    push!(Fn,ft1[1])
                    push!(Fn,ft1[2])   
                end
                if output_type == "below" || output_type == "full"
                    ft2 = TriangleFace{Int64}(f[indNext2],indMap[ind_e[indNext2]],indMap[ind_e[indNext1]])
                    push!(Fn,ft2)
                end
            end
        end
    elseif all(LV) 
        if output_type == "below" || output_type == "full"
            Fn = deepcopy(F) # Vector{QuadFace{Int64}}(,1)
            Vn = deepcopy(V) # Vector{Point3{Float64}}()
        else
            Fn = [QuadFace{Int64}] # Vector{QuadFace{Int64}}(,1)
            Vn = [Point3{Float64}] # Vector{Point3{Float64}}()
        end
    else ~any(LV) # all above
        if output_type == "above" || output_type == "full"
            Fn = deepcopy(F) # Vector{QuadFace{Int64}}(,1)
            Vn = deepcopy(V) # Vector{Point3{Float64}}()
        else
            Fn = [QuadFace{Int64}] # Vector{QuadFace{Int64}}(,1)
            Vn = [Point3{Float64}] # Vector{Point3{Float64}}()
        end
    end
    # Fn,Vn,_ = remove_unused_vertices(Fn,Vn)
    return Fn,Vn
end


function trisurfslice2(F,V,n = (0.0,0.0,1.0), p = mean(V,dims=1); snapTolerance = 0, output_type="full")

    # Compute dot product with slicing vector/plane
    d = map(v-> dot(n,v.-p),V)
    if snapTolerance != 0.0
        d[abs.(d).<snapTolerance] .= 0
    end
    LV = d.<0 # Boolean for points under the plane

    if !all(LV) && any(LV)
        if output_type == "above"
            indKeep_F = findall([~any(LV[f]) for f ∈ F])
        elseif output_type == "below"
            indKeep_F = findall([all(LV[f]) for f ∈ F])
        elseif output_type == "full"
            indKeep_F = findall([all(LV[f]) || ~any(LV[f]) for f ∈ F])
        end

        E = meshedges(F)
        E_uni,indUni,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        con_E2F = con_edge_face(F,E_uni) # EDGE-FACE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse) # FACE-EDGE connectivity          

        # Compute edge intersection points
        Vp = Vector{Point{3,Float64}}()
        ind_cut_edge = Vector{Int64}()        
        LF_cut = fill(false,length(F)) # Bool for faces to keep
        below_count_E_uni = Vector{Int64}(undef,length(E_uni))        
        for q ∈ eachindex(E_uni) # Loop over all unique edges
            e = E_uni[q] # The current edge
            below_count_E_uni[q] = sum(LV[e])
            if below_count_E_uni[q]==1 # If only 1 point is under
                pe = V[e[1]] # Vector origin
                ne = V[e[2]].-pe # Edge vector
                de = dot(n,ne) # slope in n direction
                f = d[e[1]]/de 
                push!(Vp,pe .- f.*ne) # Add new edge intersection point
                push!(ind_cut_edge,q) # Add edge index (referring to non-unique set)                
                LF_cut[con_E2F[q]] .= true
            end    
        end

        # Formulate edge to intersection point index vector
        indMap = zeros(Int64,length(E_uni))
        indMap[ind_cut_edge] = length(V)+1:1:length(ind_cut_edge)+length(V)        
        
        # Append new points to V    
        Vn = append!(deepcopy(V),Vp)         
        Fn = Vector{eltype(F)}() # F[LF_keep] # Faces to keep   
        for q ∈ eachindex(F)
            f = F[q]
            if LF_cut[q]                
                ind_e = con_F2E[q] # Edges for current face 
                below_counts = below_count_E_uni[ind_e]
                
                ind2 = findfirst(below_counts.==2)
                if ~isnothing(ind2) # Found two under                    
                    indNext1 = wrapindex(ind2+1,3)
                    indNext2 = wrapindex(ind2+2,3)                    
                    if output_type == "below" || output_type == "full"
                        fq = QuadFace{Int64}(f[ind2],f[indNext1],indMap[ind_e[indNext1]],indMap[ind_e[indNext2]])
                        ft1 = quad2tri([fq],Vn; convert_method = "angle")            
                        push!(Fn,ft1[1])
                        push!(Fn,ft1[2])
                    end
                    if output_type == "above" || output_type == "full"
                        ft2 = TriangleFace{Int64}(f[indNext2],indMap[ind_e[indNext2]],indMap[ind_e[indNext1]])
                        push!(Fn,ft2)
                    end
                else # Found two above            
                    ind2 = findfirst(below_counts.==0)
                    indNext1 = wrapindex(ind2+1,3)
                    indNext2 = wrapindex(ind2+2,3)
                    if output_type == "above" || output_type == "full"
                        fq = QuadFace{Int64}(f[ind2],f[indNext1],indMap[ind_e[indNext1]],indMap[ind_e[indNext2]])
                        ft1 = quad2tri([fq],Vn; convert_method = "angle")  
                        push!(Fn,ft1[1])
                        push!(Fn,ft1[2])   
                    end
                    if output_type == "below" || output_type == "full"
                        ft2 = TriangleFace{Int64}(f[indNext2],indMap[ind_e[indNext2]],indMap[ind_e[indNext1]])
                        push!(Fn,ft2)
                    end
                end                
            else
                if output_type == "full"            
                    push!(Fn,f)
                elseif output_type == "below"            
                    if all(LV[f])
                        push!(Fn,f)
                    end
                elseif output_type == "above"            
                    if ~any(LV[f])
                        push!(Fn,f)
                    end
                end
            end
        end
    end
    # elseif all(LV) 
    #     if output_type == "below" || output_type == "full"
    #         Fn = deepcopy(F) # Vector{QuadFace{Int64}}(,1)
    #         Vn = deepcopy(V) # Vector{Point3{Float64}}()
    #     else
    #         Fn = [QuadFace{Int64}] # Vector{QuadFace{Int64}}(,1)
    #         Vn = [Point3{Float64}] # Vector{Point3{Float64}}()
    #     end
    # else ~any(LV) # all above
    #     if output_type == "above" || output_type == "full"
    #         Fn = deepcopy(F) # Vector{QuadFace{Int64}}(,1)
    #         Vn = deepcopy(V) # Vector{Point3{Float64}}()
    #     else
    #         Fn = [QuadFace{Int64}] # Vector{QuadFace{Int64}}(,1)
    #         Vn = [Point3{Float64}] # Vector{Point3{Float64}}()
    #     end
    # end
    # Fn,Vn,_ = remove_unused_vertices(Fn,Vn)
    return Fn,Vn
end

function trisurfslice3(F,V,n = (0.0,0.0,1.0), p = mean(V,dims=1); snapTolerance = 0, output_type="full")

    # Compute dot product with slicing vector/plane
    d = map(v-> dot(n,v.-p),V)
    if snapTolerance != 0.0
        d[abs.(d).<snapTolerance] .= 0
    end
    LV = d.<0 # Boolean for points under the plane

    if !all(LV) && any(LV)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        con_E2F = con_edge_face(F,E_uni) # EDGE-FACE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse) # FACE-EDGE connectivity          

        # Compute edge intersection points
        below_count_E_uni = map(e-> sum(LV[e]),E_uni)
        ind_cut_edge = findall(below_count_E_uni.==1)
        Vp = Vector{Point{3,Float64}}(undef,length(ind_cut_edge))        
        LF_cut = fill(false,length(F)) # Bool for faces to keep        
        for q ∈ eachindex(ind_cut_edge) # Loop over all unique edges
            e = E_uni[ind_cut_edge[q]] # The current edge                       
            pe = V[e[1]] # Vector origin
            ne = V[e[2]].-pe # Edge vector
            de = dot(n,ne) # slope in n direction
            f = d[e[1]]/de 
            Vp[q] = pe .- f.*ne # Add new edge intersection point                
            LF_cut[con_E2F[ind_cut_edge[q]]] .= true            
        end

        # Formulate edge to intersection point index vector
        indMap = zeros(Int64,length(E_uni))
        indMap[ind_cut_edge] = length(V)+1:1:length(ind_cut_edge)+length(V)        
    
        # Append new points to V    
        Vn = append!(deepcopy(V),Vp)         
        Fn = Vector{eltype(F)}() # F[LF_keep] # Faces to keep   
        for q ∈ eachindex(F)
            f = F[q]
            if LF_cut[q]                
                ind_e = con_F2E[q] # Edges for current face 
                below_counts = below_count_E_uni[ind_e]
                
                ind2 = findfirst(below_counts.==2)
                if ~isnothing(ind2) # Found two under                    
                    indNext1 = wrapindex(ind2+1,3)
                    indNext2 = wrapindex(ind2+2,3)                    
                    if output_type == "below" || output_type == "full"
                        fq = QuadFace{Int64}(f[ind2],f[indNext1],indMap[ind_e[indNext1]],indMap[ind_e[indNext2]])                        
                        ft1 = quad2tri([fq],Vn; convert_method = "angle")            
                        push!(Fn,ft1[1])
                        push!(Fn,ft1[2])
                    end
                    if output_type == "above" || output_type == "full"
                        ft2 = TriangleFace{Int64}(f[indNext2],indMap[ind_e[indNext2]],indMap[ind_e[indNext1]])
                        push!(Fn,ft2)
                    end
                else # Found two above            
                    ind2 = findfirst(below_counts.==0)
                    indNext1 = wrapindex(ind2+1,3)
                    indNext2 = wrapindex(ind2+2,3)
                    if output_type == "above" || output_type == "full"
                        fq = QuadFace{Int64}(f[ind2],f[indNext1],indMap[ind_e[indNext1]],indMap[ind_e[indNext2]])
                        ft1 = quad2tri([fq],Vn; convert_method = "angle")  
                        push!(Fn,ft1[1])
                        push!(Fn,ft1[2])   
                    end
                    if output_type == "below" || output_type == "full"
                        ft2 = TriangleFace{Int64}(f[indNext2],indMap[ind_e[indNext2]],indMap[ind_e[indNext1]])
                        push!(Fn,ft2)
                    end
                end                
            else
                if output_type == "full"            
                    push!(Fn,f)
                elseif output_type == "below"            
                    if all(LV[f])
                        push!(Fn,f)
                    end
                elseif output_type == "above"            
                    if ~any(LV[f])
                        push!(Fn,f)
                    end
                end
            end
        end
    end
    return Fn,Vn
end






# println("-------------")
# @time Fn11,Vn11 = trisurfslice(F,V,n,p; output_type="below")
# @time Fn12,Vn12 = trisurfslice(F,V,n,p; output_type="below")
# println("-------------")
# @time Fn21,Vn21 = trisurfslice2(F,V,n,p; output_type="below")
# @time Fn22,Vn22 = trisurfslice2(F,V,n,p; output_type="below")
# println("-------------")
# @time Fn31,Vn31 = trisurfslice3(F,V,n,p; output_type="below")
# @time Fn32,Vn32 = trisurfslice3(F,V,n,p; output_type="below")
# println("-------------")
# @time Fn41,Vn41 = trisurfslice4(F,V,n,p; output_type="below")
# @time Fn42,Vn42 = trisurfslice4(F,V,n,p; output_type="below")
# println("-------------")

Fn,Vn = trisurfslice(F,V,n,p; output_type="below")
Mn = GeometryBasics.Mesh(Vn,Fn)


## Visualization
s = 1.25*maximum([maximum(map(v-> v[i],V)) - minimum(map(v-> v[i],V)) for i ∈ 1:3])

R = rotation_between(n,[0.0,0.0,1.0])
plateDim = (s,s)
plateElem = (1,1)
FG1,VG1 = quadplate(plateDim,plateElem)
MG = GeometryBasics.Mesh(VG1,FG1)
fig = Figure(size=(800,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")

stepRange = range(-s/2,s/2,50)
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 0,linewidth=30)

Mn = lift(hSlider.value) do stepIndex       
    pp = [p[1],p[2],p[3]+stepIndex]
    Fn,Vn = trisurfslice4(F,V,n,pp; output_type="below")        
    return GeometryBasics.Mesh(Vn,Fn)
end

MG = lift(hSlider.value) do stepIndex   
    pp = [p[1],p[2],p[3]+stepIndex]
    
    VGn = [GeometryBasics.Point{3, Float64}(R'*v) for v ∈ VG1] # Rotate plane
    VGn = map(v-> v.+pp,VGn) # Offset plate    
    return GeometryBasics.Mesh(togeometrybasics_points(VGn),FG1)
end

hp1 = mesh!(ax1,GeometryBasics.Mesh(V,F),color=:white, shading = FastShading, transparency=true)
hp2 = wireframe!(ax1,MG, linewidth=5, color=:red)
hp3 = poly!(ax1,Mn, strokewidth=1,color=:white, strokecolor=:blue, shading = FastShading, transparency=false)
# hp3 = normalplot(ax1,Mn)

slidercontrol(hSlider,ax1)
fig
