using Comodo
using GeometryBasics
using GLMakie
using Rotations
using Statistics
using LinearAlgebra

testCase = 3

if testCase == 1    
    # Define guide curve
    nc = 101 # Number of points on guide curve
    P = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
    P[1 ] = GeometryBasics.Point{3, Float64}( 0.0, 0.0, 0.0)
    P[2 ] = GeometryBasics.Point{3, Float64}( 1.0, 0.0, 0.0)
    P[3 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 0.0)
    P[4 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 1.0)
    Vc = nbezier(P,nc) # Get Bezier fit points
    Vc = [vc.*10 for vc in Vc]
    Vc,Sc = evenly_sample(Vc, nc)

    # Define section curves
    np = 50 # Number of section points
    f(t) = 2.0 + 0.5.*sin(3*t)
    V1 = circlepoints(f,np; dir=:acw)
    V1,_ = evenly_sample(V1, np)
    Q = RotXYZ(0.0,0.5*π,0.0) # Define a rotation tensor using Euler angles
    V1 = [Q*v for v ∈ V1] # Rotate the coordinates

    f(t) = 2.0 + 0.5*sin(3*t)
    V2 = circlepoints(f,np; dir=:acw)
    V2,_ = evenly_sample(V2, np)
    V2 = [v2 .+ Vc[end] for v2 ∈ V2] 
elseif testCase == 2
        # Define guide curve
        nc = 75 # Number of points on guide curve
        P = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        P[1 ] = GeometryBasics.Point{3, Float64}( 0.0, 0.0, 0.0)
        P[2 ] = GeometryBasics.Point{3, Float64}( 1.0, 0.0, 0.0)
        P[3 ] = GeometryBasics.Point{3, Float64}( 2.0, 0.0, 0.0)
        P[4 ] = GeometryBasics.Point{3, Float64}( 3.0, 0.0, 0.0)
        Vc = nbezier(P,nc) # Get Bezier fit points
        Vc = [vc.*10 for vc in Vc]
        Vc,Sc = evenly_sample(Vc, nc)
    
        # Define section curves
        np = 35 # Number of section points
        f(t) = 2.0 + 0.5.*sin(3*t)
        V1 = circlepoints(f,np; dir=:acw)
        V1,_ = evenly_sample(V1, np)
        Q = RotXYZ(0.0,0.5*π,0.0) # Define a rotation tensor using Euler angles
        V1 = [(Q*v) .+ Vc[1] for v ∈ V1] # Rotate the coordinates
    
        f(t) = 2.0 + 0.5*sin(3*t)
        V2 = circlepoints(f,np; dir=:acw)
        V2,_ = evenly_sample(V2, np)
        Q = RotXYZ(0.0,0.5*π,0.0) # Define a rotation tensor using Euler angles
        V2 = [(Q*v) .+ Vc[end] for v ∈ V2] 
elseif testCase == 3
    nc = 401 # Number of points on guide curve
    r = 10.0
    a = 4*π
    Vc = [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),10.0*(t/(a/2))) for t ∈ range(0,a,nc)]

    # Define section curves
    np = 50 # Number of section points

    # Section 1
    f(t) = 1.5 + 0.5.*sin(3*t)
    V1 = circlepoints(f,np; dir=:cw)
    V1,_ = evenly_sample(V1, np)
    Q = RotXYZ(0.5*π,0.0,0.0) # Define a rotation tensor using Euler angles
    V1 = [Q*v for v ∈ V1] # Rotate the coordinates

    # Ensure section is orthogonal to guide curve
    n3_1 = Q*Vec3{Float64}(0.0,0.0,-1.0)                       
    n2_1 = Q*Vec3{Float64}(1.0,0.0,0.0)   
    n1_1 = normalizevector(cross(n3_1,n2_1))
    S11 = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    n3_1 = normalizevector(normalizevector(Vc[2]-Vc[1]))                            
    n2_1 = normalizevector(cross(normalizevector(V1[1]),n3_1))
    n1_1 = normalizevector(cross(n3_1,n2_1))
    S12 = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    R = RotMatrix3{Float64}(S12\S11)
    V1 = [R*v for v ∈ V1]
    V1= [v .+ Vc[1] for v ∈ V1] 

    # Section 2
    f(t) = 4 + 1.5*sin(5*t)
    V2 = circlepoints(f,np; dir=:cw)
    V2,_ = evenly_sample(V2, np)
    Q1 = RotXYZ(0.5*π,0.0,0.0) # Define a rotation tensor using Euler angles
    Q2 = RotXYZ(0.0,-0.25*π,0.0) # Define a rotation tensor using Euler angles
    Q = Q2*Q1
    V2 = [Q*v for v ∈ V2] # Rotate the coordinates

    # Ensure section is orthogonal to guide curve
    n3_2 = Q*Vec3{Float64}(0.0,0.0,-1.0)                       
    n2_2 = Q*Vec3{Float64}(1.0,0.0,0.0)   
    n1_2 = normalizevector(cross(n3_2,n2_2))
    S21 = RotMatrix3{Float64}(mapreduce(permutedims,vcat,[n1_2,n2_2,n3_2]))

    n3_2 = normalizevector(normalizevector(Vc[2]-Vc[1]))                            
    n2_2 = normalizevector(cross(normalizevector(V1[1]),n3_1))
    n1_2 = normalizevector(cross(n3_1,n2_1))
    S22 = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    R = RotMatrix3{Float64}(S22\S21)    
    V2 = [R*v for v ∈ V2] # Rotate the coordinates
    V2 = [v .+ Vc[end] for v ∈ V2] 

elseif testCase == 4
    nc = 1001 # Number of points on guide curve    
    t = range(0,4.0*pi,nc)
    ff = 8.0 
    rb = 6.0
    rc = rb/3   
    r = [rb + rc*sin(ff*tt) for tt in t]    
    Vc = [GeometryBasics.Point{3, Float64}(r[i]*cos(t[i]),r[i]*sin(t[i]),t[i]+rc*cos(ff*t[i])) for i in eachindex(t)]
    Vc,_ = evenly_sample(Vc, nc)

    # Define section curves
    np = 150 # Number of section points

    # Section 1
    f(t) = rc/3 + rc/5 .* sin(3*t)
    V1 = circlepoints(f,np; dir=:cw)
    V1,_ = evenly_sample(V1, np)
    Q = RotXYZ(0.5*π,0.0,0.0) # Define a rotation tensor using Euler angles
    V1 = [Q*v for v ∈ V1] # Rotate the coordinates

    # Ensure section is orthogonal to guide curve
    n3_1 = Q*Vec3{Float64}(0.0,0.0,-1.0)                       
    n2_1 = Q*Vec3{Float64}(1.0,0.0,0.0)   
    n1_1 = normalizevector(cross(n3_1,n2_1))
    S11 = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    n3_1 = normalizevector(normalizevector(Vc[2]-Vc[1]))                            
    n2_1 = normalizevector(cross(normalizevector(V1[1]),n3_1))
    n1_1 = normalizevector(cross(n3_1,n2_1))
    S12 = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    R = RotMatrix3{Float64}(S12\S11)
    V1 = [R*v for v ∈ V1]
    V1= [v .+ Vc[1] for v ∈ V1] 

    # Section 2
    f(t) = rc/4 + rc/6*sin(5*t)
    V2 = circlepoints(f,np; dir=:cw)
    V2,_ = evenly_sample(V2, np)
    Q1 = RotXYZ(0.5*π,0.0,0.0) # Define a rotation tensor using Euler angles
    Q2 = RotXYZ(0.0,-0.25*π,0.0) # Define a rotation tensor using Euler angles
    Q = Q2*Q1
    V2 = [Q*v for v ∈ V2] # Rotate the coordinates

    # Ensure section is orthogonal to guide curve
    n3_2 = Q*Vec3{Float64}(0.0,0.0,-1.0)                       
    n2_2 = Q*Vec3{Float64}(1.0,0.0,0.0)   
    n1_2 = normalizevector(cross(n3_2,n2_2))
    S21 = RotMatrix3{Float64}(mapreduce(permutedims,vcat,[n1_2,n2_2,n3_2]))

    n3_2 = normalizevector(normalizevector(Vc[2]-Vc[1]))                            
    n2_2 = normalizevector(cross(normalizevector(V1[1]),n3_1))
    n1_2 = normalizevector(cross(n3_1,n2_1))
    S22 = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    R = RotMatrix3{Float64}(S22\S21)    
    V2 = [R*v for v ∈ V2] # Rotate the coordinates
    V2 = [v .+ Vc[end] for v ∈ V2] 

end



#########

# face_type=:quad
F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist=0)
F = invert_faces(F)
C = ceil.(collect(1:length(V))./(length(V1)))

# Visualization

cMap = Makie.Reverse(:Spectral)
markersize = 6
linewidth = 2

vizCase = 2
if vizCase ==1 

    CF = round.(Int64,vertex2simplexdata(F,C))

    fig = Figure(size = (1600,1600))
    ax = Axis3(fig[1, 1],aspect = :data,title="Swept lofting")

    stepRange1 = 0:maximum(C)

    hSlider1 = Slider(fig[2, 1], range = stepRange1, startvalue = 0,linewidth=30)

    # scatter!(ax, Vc,markersize=markersize,color=:black)
    hp1 = lines!(ax, Vc,linewidth=linewidth,color=:black)

    scatter!(ax, V1,markersize=markersize,color=:blue)
    hp2 = lines!(ax, V1,linewidth=linewidth,color=:blue)

    scatter!(ax, V2,markersize=markersize,color=:red)
    hp3 = lines!(ax, V2,linewidth=linewidth,color=:red)

    # hp1 = poly!(ax, GeometryBasics.Mesh(V,F), color=C, strokecolor=:black, strokewidth=0.5,transparency=false,shading = FastShading,colormap=:Spectral)
    hp1 = mesh!(ax, GeometryBasics.Mesh(V,F), color=C,transparency=false,shading = FastShading,colormap=cMap)

    on(hSlider1.value) do stepIndex1        
        hp1[1] = GeometryBasics.Mesh(V,F[CF.<=stepIndex1])  
    end

    # fileName = comododir()*"/assets/img/sweep_loft_anim_01.mp4"
    # slider2anim(fig,hSlider1,fileName; backforth=true, duration=10)

    fig 

elseif vizCase ==2 
    fig = Figure(size = (800,800))
    ax = Axis3(fig[1, 1],aspect = :data,title="Swept lofting")

    stepRange1 = -10:1:10
    hSlider1 = Slider(fig[2, 1], range = stepRange1, startvalue = 0,linewidth=30)

    # scatter!(ax, Vc,markersize=markersize,color=:black)
    hp1 = lines!(ax, Vc,linewidth=linewidth,color=:black)

    scatter!(ax, V1,markersize=markersize,color=:blue)
    hp2 = lines!(ax, V1,linewidth=linewidth,color=:blue)

    scatter!(ax, V2,markersize=markersize,color=:red)
    hp3 = lines!(ax, V2,linewidth=linewidth,color=:red)

    # hp1 = poly!(ax, GeometryBasics.Mesh(V,F), strokecolor=:black, strokewidth=0.5,color=C,transparency=false,shading = FastShading,colormap=:Spectral)
    hp1 = mesh!(ax, GeometryBasics.Mesh(V,F), color=C,transparency=false,shading = FastShading,colormap=cMap)

    on(hSlider1.value) do stepIndex1
        F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist=stepIndex1)
        F = invert_faces(F)   
        hp1[1] = GeometryBasics.Mesh(V,F)    
    end

    # fileName = comododir()*"/assets/img/sweep_loft_anim_02.mp4"
    # slider2anim(fig,hSlider1,fileName; backforth=true, duration=3)

    fig

end

