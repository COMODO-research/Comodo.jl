using Comodo
using Comodo.GLMakie
using Comodo.LinearAlgebra
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        V = [Point{3,Float64}(0.0, 0.0, 0.0), Point{3,Float64}(1.0, 0.0, 0.0), Point{3,Float64}(0.5, 1.0, 0.0)]
        F = [TriangleFace{Int}(1,2,3)]
    elseif testCase == 2
        r = 2.0
        n = 2
        F, V = geosphere(n,r; method=:Loop)  
    elseif testCase == 3    
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)        
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)    
    end

    R, P = circumcircle(F,V)
    
    # Visualization
    depth_shift = -0.01f0

    fig = Figure(size=(1200,800))    
    ax1 = AxisGeom(fig[1,1]; title = "Triangle(s) with circumcircle(s)")
    hp1 = meshplot!(ax1,F,V)
    scatter!(ax1, P, markersize=8, color = :red, depth_shift=depth_shift)

    np = 35
    m = np + 2
    VC = Vector{Point{3,Float64}}(undef, m*length(P))
    nz = Point{3,Float64}(0.0, 0.0, 1.0)
    Vc = circlepoints(1.0,np) 
    push!(Vc,Vc[1])
    push!(Vc,Point{3,Float64}(NaN,NaN,NaN))    
    
    for (i,p) in enumerate(P)           
        n = facenormal(F[i],V)      
        Q = rotation_between(nz,n)
        VC[1+(i-1)*m:i*m] = [Q*R[i]*v+P[i] for v in Vc]
    end
    lines!(ax1, VC, linewidth=1, color = :red, depth_shift=depth_shift)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end

# save("/home/kevin/Desktop/im$testCase.jpg", fig, px_per_unit = 2)