using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

GLMakie.closeall()

for testCase = 1:8
    if testCase == 1
        r = 25.25
        F,V = geosphere(5,r)  
        # V = [GeometryBasics.Point{3, Float64}(v[1],v[2],3*v[3]) for v ∈ V]  
    elseif testCase==2
        r = 25.25
        F,V = quadsphere(3,r)    
    elseif testCase==3
        r = sqrt(3)
        F,V = cube(r)    
        # F = quad2tri(F,V; convert_method = "angle")
    elseif testCase==4   
        r = 25
        s = r/10
        nc = round(Int,(2*pi*r)/s)
        d = r*2
        Vc = circlepoints(r, nc; dir=:cw)
        num_steps = round(Int,d/s)
        num_steps = num_steps + Int(iseven(num_steps))
        F, V = extrudecurve(Vc; extent=d, n=[0.0,0.0,1.0], num_steps=num_steps, close_loop=true, face_type=:tri)
        R = rand(RotXYZ)
        V = [R* v for v in V]
    elseif testCase==5 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
        # F,V = subtri(F,V,2; method = :Loop)
    elseif testCase==6 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
    elseif testCase==7
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","obj","motherChild_5k.obj")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
    elseif testCase == 8
        nSub = 4 # Number of refinement steps of the geodesic sphere
        r = 2.5 # Sphere radius
        F,V = geosphere(nSub,r) # Creating the faces and vertices of a full sphere
        VC = simplexcenter(F,V) # Finding triangle centre coordinates
        F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces using z of central coordinates
        F,V = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed
    end
        
    M = GeometryBasics.Mesh(V,F)

    K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V; growsteps=2)

    ## Visualization
    s = pointspacingmean(F,V)
    cMap = Makie.Reverse(:Spectral)
    f = 0.1

    fig = Figure(size=(800,800))

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "1st principal curvature")
    hp1 = poly!(ax1,M, color=K1, strokecolor=:white, strokewidth=0.1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(filter(!isnan,K1))).*f.*(-1,1))
    hpn1 = dirplot(ax1,V,U1; color=:black,linewidth=2,scaleval=s,style=:through)
    # scatter!(ax1,V,color=K1,colormap=cMap,colorrange = maximum(abs.(K1)).*0.1.*(-1,1),markersize=10);
    Colorbar(fig[1, 2],hp1)

    ax1 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "2nd principal curvature")
    hp2 = poly!(ax1,M, color=K2, strokecolor=:white, strokewidth=0.1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(filter(!isnan,K2))).*f.*(-1,1))
    hpn2 = dirplot(ax1,V,U2; color=:black,linewidth=2,scaleval=s,style=:through)
    # scatter!(ax1,V,color=K2,colormap=cMap,colorrange = maximum(abs.(K2)).*0.1.*(-1,1),markersize=10);
    Colorbar(fig[1, 4],hp2)

    ax1 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Mean curvature")
    hp3 = poly!(ax1,M, color=H, strokecolor=:white, strokewidth=0.1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(filter(!isnan,H))).*f.*(-1,1))
    Colorbar(fig[2, 2],hp3)

    ax1 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Gaussian curvature")
    hp4 = poly!(ax1,M, color=G, strokecolor=:white, strokewidth=0.1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(filter(!isnan,G))).*f.*(-1,1))
    Colorbar(fig[2, 4],hp4)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end