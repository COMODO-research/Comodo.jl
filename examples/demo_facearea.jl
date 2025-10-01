using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using FileIO

GLMakie.closeall()

for testCase = 1:6
    # Example geometry
    if testCase == 1
        r = 1
        F,V = geosphere(3,r)  
        # V = [GeometryBasics.Point{3, Float64}(v[1],v[2],3*v[3]) for v ∈ V]  
        # V = [v + Point{3,Float64}(2.0, 2.0, 2.0,) for v in V]

        println("Theoretical area: " * string(4*pi*r^2))    
    elseif testCase == 2
        r = 1.0
        F,V = quadsphere(3,r)    
        println("Theoretical area: " * string(4*pi*r^2))    
    elseif testCase == 3
        r = 2*sqrt(3)/2
        F,V = cube(r)        
        println("Theoretical area: " * string(6*(2*(r./(sqrt(3))))^2))
    elseif testCase == 4 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)    
    elseif testCase == 5 # Merged STL for single object
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)
    elseif testCase == 6
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","obj","motherChild_5k.obj")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
    end

    M = GeometryBasics.Mesh(V,F)

    A = facearea(F,V)
    println("Total mesh area: " *string(sum(A)))

    Fn,Vn = separate_vertices(F,V)

    ## Visualization
    fig = Figure(size=(800,800))
    ax1 = AxisGeom(fig[1, 1], title = "Face area")
    hp1 = meshplot!(ax1, Fn, Vn, color=simplex2vertexdata(Fn,A), strokewidth=0.5)
    Colorbar(fig[1, 2],hp1)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end