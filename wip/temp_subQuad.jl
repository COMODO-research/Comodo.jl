using Gibbon
using GLMakie
using GeometryBasics

## Define example input
r = 0.5 #radius
M = platonicsolid(2,r) # Get an example quadrilateral mesh (for a cube in this case)
V = coordinates(M)
F = faces(M)

# function subQuadt(F,V,n)
#     if n==0
#         return F,V
#     elseif n==1

#         # Get edges
#         E = meshEdges(F) # Non-unique edges

#         Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
#         con_F2E = con_face_edge(F,Eu,indReverse)

#         # Define vertices
#         Ve = midPoints(Eu,V) # Mid edge points
#         Vf = midPoints(F,V)  # Mid face points
#         Vn = [V;Ve;Vf] # Joined point set

#         # Define faces
#         Fn = Vector{QuadFace{Int64}}(undef,length(F)*4)        
#         nv = length(V)
#         ne = length(Eu)
#         nf = length(F)
#         for q ∈ eachindex(F)
#             i = 1 + (q-1)*4
#             for ii = 0:1:3                
#                 Fn[i+ii] = QuadFace{Int64}([F[q][ii+1], con_F2E[q][ii+1]+nv, q+nv+ne, con_F2E[q][1+mod(3+ii,4)]+nv])
#             end            
#         end

#         return Fn,Vn

#     elseif n>1
#         for _ =1:1:n
#             F,V = subQuadt(F,V,1)
#         end

#         return F,V
#     end
# end

## Refine mesh using `subQuad` and the default "linear" method
Fn1,Vn1=subQuad(F,V,1) # Split once 

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=1")
wireframe!(ax1,M, linewidth=8,color=:red, overdraw=false)
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=3,color=:white,shading=FastShading,transparency=false)


Legend(fig[1, 2],[hp1,hp2],["Initial","Refined"])

fig