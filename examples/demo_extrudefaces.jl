using Comodo
using GLMakie
using GeometryBasics
using Rotations
using Statistics
using LinearAlgebra

testCase = 1

if testCase ==1 
    plateDim = [10.0,10.0]
    plateElem = [5,5]
    F,V = quadplate(plateDim,plateElem)

    t = 6
    n = 4
    direction=:both
elseif testCase == 2
    plateDim = [10.0,10.0]
    plateElem = [5,5]
    F,V = quadplate(plateDim,plateElem)
    F = quad2tri(F,V)

    t = 6
    n = 4
    direction=:positive
elseif testCase == 3
    nc = 101 # Number of points on guide curve
    r = 8
    a = 3.0*π
    Vc = [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),12.0*(t/(a/2))) for t ∈ range(0,a,nc)]

    # Define section curves
    np = 50 # Number of section points

    # Section 1
    f(t) = 3 + 0.25.*sin(3*t)
    V1 = circlepoints(f,np; dir=:cw)
    V1,_ = evenly_sample(V1, np)
    V1_ori = deepcopy(V1)

    # Ensure section is orthogonal to guide curve
    n3 = Vec3{Float64}(0.0,0.0,1.0)                       
    n2 = Vec3{Float64}(0.0,1.0,0.0)   
    n1 = normalizevector(cross(n2,n3))
    S11 = mapreduce(permutedims,vcat,[n1,n2,n3])

    n3 = normalizevector(normalizevector(Vc[2]-Vc[1]))                            
    n1 = normalizevector(V1[1])
    n2 = normalizevector(cross(n1,n3))
    n1 = normalizevector(cross(n3,n2))
    S12 = mapreduce(permutedims,vcat,[n1,n2,n3])

    R = RotMatrix3{Float64}(S12\S11)
    V1 = [R*v for v ∈ V1]
    V1= [v .+ Vc[1] for v ∈ V1] 

    # Section 2
    f(t) = 5 + 0.5*sin(5*t)
    V2 = circlepoints(f,np; dir=:cw)
    V2,_ = evenly_sample(V2, np)    
    V2_ori = deepcopy(V2)

    # Q1 = RotXYZ(0.5*π,0.0,0.0) # Define a rotation tensor using Euler angles
    # Q2 = RotXYZ(0.0,-0.25*π,0.0) # Define a rotation tensor using Euler angles
    # Q = Q2*Q1
    # V2 = [Q*v for v ∈ V2] # Rotate the coordinates

    # Ensure section is orthogonal to guide curve
    n3 = Vec3{Float64}(0.0,0.0,1.0)                       
    n2 = Vec3{Float64}(0.0,1.0,0.0)   
    n1 = normalizevector(cross(n2,n3))
    S21 = RotMatrix3{Float64}(mapreduce(permutedims,vcat,[n1,n2,n3]))

    n3 = normalizevector(normalizevector(Vc[end]-Vc[end-1]))                            
    n1 = normalizevector(Point{3,Float64}(Vc[end][1],Vc[end][2],0))
    n2 = normalizevector(cross(n1,n3))
    n1 = normalizevector(cross(n3,n2))
    S22 = mapreduce(permutedims,vcat,[n1,n2,n3])

    R = RotMatrix3{Float64}(S22\S21)    
    V2 = [R*v for v ∈ V2] # Rotate the coordinates
    V2 = [v .+ Vc[end] for v ∈ V2] 

#########

    # face_type=:quad
    F,V = sweeploft(Vc,V1,V2; face_type=:tri, num_twist=0)

    t = 1
    n = 5
    direction=:positive
end



N = vertexnormal(F,V; weighting=:area)
E, Ve = extrudefaces(F,V; extent=t, direction=direction, num_steps=n)

FE = element2faces(E)

## Visualize mesh

strokewidth=1

fig = Figure(size=(1200,1200))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Input surface mesh")


hp2=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=strokewidth,color=:white, shading = FastShading,transparency=false)
# normalplot(ax1,F,V; scaleval=0.25)

# hp2 = scatter!(ax1, Ve,markersize=15,color=:orange)
# normalplot(ax1,FE,Ve; scaleval=0.25)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Extruded mesh")

if isa(FE,Tuple)
    FE1s,V1s = separate_vertices(FE[1],Ve)
    FE2s,V2s = separate_vertices(FE[2],Ve)

    ind = reduce(vcat,FE[1])
    C1 = repeat(collect(1:n),inner=length(V))
    Cs1 = C1[ind]

    ind = reduce(vcat,FE[2])
    C2 = repeat(collect(1:n),inner=length(V))
    Cs2 = C2[ind]

    hp2=poly!(ax2,GeometryBasics.Mesh(V1s,FE1s), strokewidth=strokewidth,color=Cs1, shading = FastShading,transparency=false,colormap=:Spectral)    
    hp3=poly!(ax2,GeometryBasics.Mesh(V2s,FE2s), strokewidth=strokewidth,color=Cs2, shading = FastShading,transparency=false,colormap=:Spectral)
#     normalplot(ax2,FE[1],Ve; scaleval=0.25)
#     normalplot(ax2,FE[2],Ve; scaleval=0.25)
else
    ind = reduce(vcat,FE)
    C = repeat(collect(1:n),inner=length(V))
    Cs = C[ind]
    FEs,Vs = separate_vertices(FE,Ve)
    hp2=poly!(ax2,GeometryBasics.Mesh(Vs,FEs), strokewidth=strokewidth,color=Cs, shading = FastShading,transparency=true,colormap=:Spectral)
    normalplot(ax2,FE,Ve; scaleval=0.25)

    poly!(ax2,GeometryBasics.Mesh(V,F), strokewidth=strokewidth,color=:lightgreen, shading = FastShading,transparency=false)
end

fig

# save(comododir()*"/assets/img/extrudeface.png",fig,px_per_unit = 4)
