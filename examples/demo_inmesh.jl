using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics
using FileIO

GLMakie.closeall()

for testCase in 1:5
    if testCase == 1
        F, V = geosphere(2, 15.0)
        nSteps = 10
        w = 16
        xr,yr,zr = ntuple(_->range(0.0, w,nSteps),3)
    elseif testCase == 2
        r1 = 10.0
        r2 = 5.0
        F1, V1 = geosphere(2, r1) # Big sphere 
        F2, V2 = geosphere(2, r2) # Small inner sphere 
        invert_faces!(F2) # Invert inner 
        F,V = joingeom(F1, V1, F2, V2)

        nSteps = 15
        w = r1+r1/10.0
        xr,yr,zr = ntuple(_->range(0.0, w,nSteps),3)
    elseif testCase == 3
        r = 1.0
        nc = 16
        t = range(2.0*π-(2.0*π/nc),0,nc)
        Vc = [Point3{Float64}(2.0+cos(tt),0.0,sin(tt)) for tt ∈ t]
        n = Vec{3, Float64}(0.0,0.0,1.0)
        num_steps = 24
        close_loop = true

        F,V = revolvecurve(Vc; extent = (2*pi - (2*pi/num_steps)), direction = :negative, 
        n = n, num_steps = num_steps, periodicity = (true,true), face_type = :tri)

        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (r/4.0, r/4.0, r/4.0)
        xr = 0.0:voxelSize[2]:maxV[1]
        yr = 0.0:voxelSize[1]:maxV[2]
        zr = minV[3]:voxelSize[3]:maxV[3] 
    elseif testCase == 4
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)

        pointSpacing = pointspacingmean(F,V)
        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (pointSpacing, pointSpacing, pointSpacing)
        xr = 0.0:voxelSize[2]:maxV[1]
        yr = 0.0:voxelSize[1]:maxV[2]
        zr = minV[3]:voxelSize[3]:maxV[3] 
    elseif testCase == 5
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
        V .-= mean(V) #Point{3, Float64}(0.0, 10.0, 20.0)

        F2, V2= geosphere(2, 20.0)
        # V2 .-= mean(V) #Point{3, Float64}(0.0, 10.0, 20.0)
        invert_faces!(F2)
        F,V = joingeom(F,V,F2,V2)

        pointSpacing = pointspacingmean(F,V)
        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (pointSpacing, pointSpacing, pointSpacing)
        xr = 0.0:voxelSize[2]:maxV[1]
        yr = 0.0:voxelSize[1]:maxV[2]
        zr = minV[3]:voxelSize[3]:maxV[3]  
    end
    
    PG = collect(gridpoints(xr, yr, zr))
    B = inmesh(F, V, PG)

    fig = Figure(size=(1200,800))   
    ax1 = AxisGeom(fig[1,1])
    hp1 = meshplot!(ax1, F, V, color=(:white, 0.5), strokewidth=0.0, transparency=true)
    hp2 = scatter!(ax1, PG[B], markersize=10, color=:green)
    hp2 = scatter!(ax1, PG[.!B], markersize=5, color=:red)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end