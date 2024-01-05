using GLMakie
using GeometryBasics

function hexFaces(r=1.0)

    # Create vertices       
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,8)
    s = r/sqrt(3.0)
    V[1 ] = GeometryBasics.Point{3, Float64}( -s, -s, -s)
    V[2 ] = GeometryBasics.Point{3, Float64}( -s,  s, -s)
    V[3 ] = GeometryBasics.Point{3, Float64}(  s,  s, -s)
    V[4 ] = GeometryBasics.Point{3, Float64}(  s, -s, -s)
    V[5 ] = GeometryBasics.Point{3, Float64}( -s, -s,  s)
    V[6 ] = GeometryBasics.Point{3, Float64}( -s,  s,  s)
    V[7 ] = GeometryBasics.Point{3, Float64}(  s,  s,  s)
    V[8 ] = GeometryBasics.Point{3, Float64}(  s, -s,  s)
        
    # Create faces
    F = Vector{QuadFace{Int64}}(undef,6)
    F[1 ] = QuadFace{Int64}(1,2,3,4)
    F[2 ] = QuadFace{Int64}(8,7,6,5)
    F[3 ] = QuadFace{Int64}(5,6,2,1)
    F[4 ] = QuadFace{Int64}(6,7,3,2)    
    F[5 ] = QuadFace{Int64}(7,8,4,3)    
    F[6 ] = QuadFace{Int64}(8,5,1,4)    

    return GeometryBasics.Mesh(V,F), V, F
end

M,V,F = hexFaces(1.0)  
CF_type =[1,2,3,4,5,6] # Face colors

#Visualize mesh
fig = Figure()

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Single color")
hp1=poly!(ax1,M,strokewidth=5,color=:red, transparency=false, overdraw=false)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Per face color attempt 1")
hp2=poly!(ax2,[GeometryBasics.Mesh(V,f) for f ∈ eachrow(F)],strokewidth=5,color=CF_type, transparency=true, overdraw=false,colormap = Reverse(:Spectral), colorrange=(0,6))

ax3=Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Per face color attempt 2")
hp3=[]
for i ∈ eachindex(F)    
    hp3=poly!(ax3,GeometryBasics.Mesh(V,F[i]),strokewidth=5,color=CF_type[i], transparency=true, overdraw=false,colormap = Reverse(:Spectral), colorrange=(0,6))
end

fig
