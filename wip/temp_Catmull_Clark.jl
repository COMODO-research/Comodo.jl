using Gibbon
using GLMakie
using GeometryBasics

r=1.0 #radius
M=platonicsolid(2,r)
V = coordinates(M)
F = faces(M)

function subQuad(F,V,n)
    if n==0
        return F,V
    elseif n==1
        E = meshEdges(F) # Non-unique edges
        Eu, ~, indE_uni = unique_simplices(E,V) # Unique edges
        
        Ve = midPoints(Eu,V) # Mid edge points
        Vf = midPoints(F,V)  # Mid face points
        Vn = [V;Ve;Vf] # Joined point set

        # Define faces
        Fn = Vector{QuadFace{Int64}}(undef,length(F)*4)        
        nv = length(V)
        ne = length(Eu)
        for q ∈ eachindex(F)
            i = 1 + (q-1)*4
            for ii = 0:1:3
                Fn[i+ii] = QuadFace{Int64}([F[q][ii+1],indE_uni[i+ii]+nv,nv+ne+q,indE_uni[i+mod(3+ii,4)]+nv])                
            end            
        end
        return Fn,Vn
    elseif n>1
        for _ =1:1:n
            F,V = subQuad(F,V,1)
        end

        return F,V
    end
end

n = 2

Fn,Vn=subQuad(F,V,3)

#Visualize mesh
fig = Figure(size=(1600,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined quadrilateral mesh")
hp1=wireframe!(ax1,M, linewidth=8,color=:red, overdraw=false)
hp2=poly!(ax1,GeometryBasics.Mesh(Vn,Fn), strokewidth=3,color=:white,shading=FastShading,transparency=false)
Legend(fig[1, 2],[hp1,hp2],["Initial","Refined"])

fig