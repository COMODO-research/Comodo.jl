using Comodo
using Comodo.GLMakie
using Comodo.LinearAlgebra
using Comodo.Statistics
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

GLMakie.closeall()
for testCase = 1:3
    if testCase == 1
        F, V = geosphere(2, 1.0)
        V_mean = Point{3,Float64}(5.0, 2.0, 3.0)
        V = [Point{3,Float64}(4.0*v[1], 2.0*v[2], 0.75*v[3]) for v in V]

        # Define a rotation tensor using Euler angles
        Q = RotXYZ(-0.25*π, 0.25*π, 0.25*pi)

        # Rotate the coordinates
        V = [Point{3, Float64}(Q*v) for v ∈ V] 
        V .+= V_mean

        S = pointsvd(V)    
        R = S.V

    elseif testCase == 2
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","femur.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)  
        Q = RotXYZ(0.0*pi, 0.0*pi, 0.25*pi)
        V = [Point{3, Float64}(Q*v) for v ∈ V] 

        F, V = mergevertices(F, V)

        S = pointsvd(V)    
        R = S.V
    elseif testCase == 3
        
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)          
        F, V = mergevertices(F, V)

        for n = 1:3
            B = [mean(V[f])[1]> 10.0 for f in F]               
            F, V = subtri_dual_local(F, V, B; smooth=true)
        end

        # Face area weighted
        VF = simplexcenter(F, V)
        A = facearea(F,V)
        S = pointsvd((A/sum(A)).*VF)    
        R = S.V
    end 

    pc = surface_centroid(F, V)

    #Rotate using svd directions
    Vr = [Point{3, Float64}(S.Vt*v) for v ∈ V] 
    
    # Visualize mesh
    s = norm(maxp(Vr)-minp(Vr)) # Vector size
    fig = Figure(size = (1600,800))

    ax1 = AxisGeom(fig[1, 1], title = "SVD derived orientations")
    hp1 = meshplot!(ax1, F, V; color=(:white, 0.5), strokewidth=0.5, transparency=true)
    arrows3d!(ax1, pc, Point{3,Float64}(R[:,1]), color = :red, lengthscale = s/2)
    arrows3d!(ax1, pc, Point{3,Float64}(R[:,2]), color = :green, lengthscale = s/2)
    arrows3d!(ax1, pc, Point{3,Float64}(R[:,3]), color = :blue, lengthscale = s/2)
    
    # normalplot(ax1, F,V )

    ax2 = AxisGeom(fig[1, 2], title = "Rotated")
    hp1 = meshplot!(ax2, F, Vr; strokewidth=0.5)
    # normalplot(ax2, F,Vr)

    screen = display(GLMakie.Screen(), fig)
end

