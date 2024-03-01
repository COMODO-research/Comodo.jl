using Comodo
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra

# Example curve
r = 1
nc = 16
Vc = circlepoints(r,nc;dir=:cw)

d = 3.0
n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))
pointSpacing =0.25 
s = 1

function extrudecurve(V1,d; s=1, n=Point{3, Float64}(0.0,0.0,1.0),num_steps=missing,close_loop=false,face_type="quad")
    if ismissing(num_steps)
        num_steps = ceil(Int64,d/pointspacingmean(V1))
        if face_type==:tri
            num_steps = num_steps + Int64(iseven(num_steps)) # Force uneven
        end
    end
    if s==1 # Allong n from V1
        p = d.*n
    elseif s==-1 # Against n from V1
        p = -d.*n
    elseif s==0 # Extrude both ways from V1
        p = d.*n
        V1 = [(eltype(V1))(v.-p./2) for v ∈ V1] #Shift V1 in negative direction
    end
    V2 = [(eltype(V1))(v.+p) for v ∈ V1]  
    return loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type=face_type)
end

#   num_loft = ceil(Int64,d/pointSpacing)
F1,V1 = extrudecurve(Vc,d;s=1, close_loop=true,face_type=:quad)
F2,V2 = extrudecurve(Vc,d;s=0, close_loop=true,face_type=:tri)
F3,V3 = extrudecurve(Vc,d;s=-1, close_loop=true,face_type=:tri_slash)

M1 = GeometryBasics.Mesh(V1,F1)
M2 = GeometryBasics.Mesh(V2,F2)
M3 = GeometryBasics.Mesh(V3,F3)

## Visualization
fig = Figure(size=(1200,1200))

ax1 = Axis3(fig[1, 1], aspect = :data, limits=(-r,r,-r,r,-d,d),xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Extruded s=1, face_type="quad" """)
hp1 = lines!(ax1,Vc,color=:red,linewidth=6, transparency=true, depth_shift=-1.0f-3)
hp2 = poly!(ax1,M1, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, limits=(-r,r,-r,r,-d,d), xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Extruded s=0, face_type="tri" """)
hp1 = lines!(ax2,Vc,color=:red,linewidth=6, transparency=true, depth_shift=-1.0f-3)
hp2 = poly!(ax2,M2, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, limits=(-r,r,-r,r,-d,d), xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Extruded s=-1, face_type="tri_slash" """)
hp1 = lines!(ax3,Vc,color=:red,linewidth=6, transparency=true, depth_shift=-1.0f-3)
hp2 = poly!(ax3,M3, strokewidth=1,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

# Legend(fig[1, 4],[hp1,hp2],["Input curve","Surface"])

fig
