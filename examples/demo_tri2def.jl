using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Rotations
using Comodo.Statistics

#=
This demo shows the use of the `tri2def` function to compute the deformation 
gradient tensor for a triangle or set of triangles.
=#

GLMakie.closeall()

for testCase = 1:3
    # Example initial mesh 
    if testCase == 1 
        P1 = [Point{3,Float64}(0.0, 0.0, 0.0), 
            Point{3,Float64}(2.0, 0.0, 0.0), 
            Point{3,Float64}(2.0, 2.0, 0.0),
            Point{3,Float64}(0.0, 2.0, 0.0),]
        TRI = [TriangleFace{Int}(1,2,3),TriangleFace{Int}(3,4,1)]
    elseif testCase == 2
        TRI, P1 = geosphere(2,5.0)
    elseif testCase == 3 
        pointSpacing = 0.5
        boxDim = [2.0, 2.0, 2.0] # Dimensions for the box in each direction
        TRI, P1, _ = tribox(boxDim,pointSpacing)
    end
    # Define true global deformation gradient tensor

    # Stretches
    λ₁ = 2.0
    λ₂ = 0.5
    λ₃ = 1.5

    # Right stretch tensor
    U = [λ₁  0.0 0.0;
        0.0 λ₂  0.0;
        0.0 0.0  λ₃]

    # Rotation tensor
    Q = RotXYZ(0.25*pi,0.0*pi, 0.25*π)

    # Deformation gradient tensor
    F = Q*U

    # Create deformed coordinate set
    P2 = [Point{3, Float64}(F*p) for p in P1]

    # Compute triangle deformation gradient tensors 
    F_set = tri2def(TRI, P1, P2)

    # -----------------------------------------------------------------------------
    # Visualization
    N1 = Vector{Point{3,Float64}}(undef,length(TRI))
    N2 = Vector{Point{3,Float64}}(undef,length(TRI))
    N3 = Vector{Point{3,Float64}}(undef,length(TRI))
    M1 = Vector{Point{3,Float64}}(undef,length(TRI))
    M2 = Vector{Point{3,Float64}}(undef,length(TRI))
    M3 = Vector{Point{3,Float64}}(undef,length(TRI))
    L1 = Vector{Float64}(undef,length(TRI))
    L2 = Vector{Float64}(undef,length(TRI))
    L3 = Vector{Float64}(undef,length(TRI))
    for (i,F) in enumerate(F_set)
        Un, Vn, Qn, W, Σ, R = polarDecomposition(F)

        # Current normal direction 
        n = facenormal(TRI[i],P1) 
        for c in 1:3
            d = dot(R[:,c], n)         
            if d>0.5 # Allignes with normal
                break
            elseif d<-0.5 # Ati-allign with normal
                W .*= -1.0
                R .*= -1.0
                break
            end
        end    
        n1 = Point{3,Float64}(R[:,1])
        n2 = Point{3,Float64}(R[:,2])
        n3 = Point{3,Float64}(R[:,3])

        m1 = Point{3,Float64}(W[:,1])
        m2 = Point{3,Float64}(W[:,2])
        m3 = Point{3,Float64}(W[:,3])

        N1[i] = n1
        N2[i] = n2
        N3[i] = n3
        M1[i] = m1
        M2[i] = m2
        M3[i] = m3
        L1[i] = Σ[1]
        L2[i] = Σ[2]
        L3[i] = Σ[3]
    end

    P1c = simplexcenter(TRI,P1)
    P2c = simplexcenter(TRI,P2)

    TRI_s, P2s = separate_vertices(TRI, P2)
    L1s = simplex2vertexdata(TRI_s, L1)
    L2s = simplex2vertexdata(TRI_s, L2)
    L3s = simplex2vertexdata(TRI_s, L3)
    cRange = (minimum(diag(U)), maximum(diag(U)))
    cMap = :Spectral
    veclength = pointspacingmean(TRI,P1)/2.0

    fig = Figure(size=(1200,800))

    ax1 = AxisGeom(fig[1, 1], title = "Initial")
    hp1 = meshplot!(ax1, TRI, P1)
    dirplot(ax1, P1c, N1; scaleval=veclength, linewidth=1.0, color=:red, style=:through)
    dirplot(ax1, P1c, N2; scaleval=veclength, linewidth=1.0, color=:green, style=:through)
    dirplot(ax1, P1c, N3; scaleval=veclength, linewidth=1.0, color=:blue, style=:through)


    ax2 = AxisGeom(fig[1, 2], title = "Deformed λ₁")
    hp2 = meshplot!(ax2, TRI_s, P2s, strokecolor=:white, color=L1s, colorrange=cRange, colormap=cMap)
    dirplot(ax2, P2c, M1; scaleval=veclength, linewidth=1.0, color=:black, style=:through)

    ax3 = AxisGeom(fig[2, 1], title = "Deformed λ₂")
    hp3 = meshplot!(ax3, TRI_s, P2s, strokecolor=:white, color=L2s, colorrange=cRange, colormap=cMap)
    dirplot(ax3, P2c, M2; scaleval=veclength, linewidth=1.0, color=:black, style=:through)

    ax4 = AxisGeom(fig[2, 2], title = "Deformed λ₃")
    hp4 = meshplot!(ax4, TRI_s, P2s, strokecolor=:white, color=L3s, colorrange=cRange, colormap=cMap)
    dirplot(ax4, P2c, M3; scaleval=veclength, linewidth=1.0, color=:black, style=:through)

    Colorbar(fig[:,3], hp2)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end