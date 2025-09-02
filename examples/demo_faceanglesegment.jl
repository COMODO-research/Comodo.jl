using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1 
        # A "zig-zag" curve with known angles for testing
        r = 2.0
        B = [170.0,150.0,135.0,90.0,75.0,60.0,45.0,-45.0,-60.0,-75,-90.0,-135.0,-150.0,-170.0]
        V = Vector{Point{3, Float64}}()
        push!(V,[0.0, 0.0, 0.0])
        for b in B
            push!(V,V[end]+[r, 0.0, 0.0])
            push!(V,V[end]+[r*cosd(b),r*sind(b), 0.0])
        end
        push!(V,V[end]+[r, 0.0, 0.0])

        F,V = extrudecurve(V; extent=1, direction=:positive, n=Vec{3, Float64}(0.0,0.0,-1.0), close_loop=false, face_type=:quad)
    elseif testCase == 2 
        # An imported STL based geometry of an engineering part with various flat faces and some known angles (e.g. 0, 45, and 90 degrees)
        # This model is relatively evenly sampled and of a relatively low resolution
        fileName_mesh = joinpath(comododir(),"assets","stl","spur_gear_01.stl")
        M = load(fileName_mesh)
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_,_ = mergevertices(F,V)
    elseif testCase == 3 # Same as 2 but higher angular/general resolution
        # An imported STL based geometry of an engineering part with various flat faces and some known angles (e.g. 0, 45, and 90 degrees)
        # This model features a non-homogeneous mesh (includes sharp triangles) and is of a relatively high resolution
        fileName_mesh = joinpath(comododir(),"assets","stl","spur_gear_02.stl")
        M = load(fileName_mesh)
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_,_ = mergevertices(F,V)
    end

    A,E,con_E2F = edgefaceangles(F,V; deg=true)

    G = faceanglesegment(F,V; deg=true, angleThreshold = 22.5, indStart = 1)

    ## Visualization
    linewidth = 5
    c = cgrad(:Spectral,maximum(G),categorical = true)
    A,E,con_E2F = edgefaceangles(F,V; deg=true)

    Lp = .~isnan.(A) .&& abs.(A).>30
    Ap = A[Lp]
    Ep = E[Lp]

    En,Vn_E = separate_vertices(Ep,V)
    An = simplex2vertexdata(En,Ap,Vn_E)


    Fn,Vn = separate_vertices(F,V)
    Gn = simplex2vertexdata(Fn,G,Vn)## Visualization

    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Imported mesh")
    hp1 = meshplot!(ax1, Fn, Vn, strokewidth=0.1)
    # normalplot(ax1,F,V)
    # hp_A = wireframe!(ax1,GeometryBasics.Mesh(Vn_E,En),linewidth=linewidth, transparency=false, color=An,colormap=:Spectral,colorrange = (-120, 120))
    # Colorbar(fig[1, 2],hp_A, label = "Angles")

    ax2 = AxisGeom(fig[1, 2], title = "Groups")
    hp2 = meshplot!(ax2, Fn, Vn, color=Gn, colormap=c)#,strokecolor=:black,strokewidth=0.25)
    # normalplot(ax1,F,V)

    hp3 = edgeplot!(ax2, En, Vn_E, linewidth=linewidth)

    Colorbar(fig[1, 3],hp2, label = "Groups")

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end