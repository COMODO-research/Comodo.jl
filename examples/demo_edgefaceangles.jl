using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

testCase = 1

if testCase==1 
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

    F,V = extrudecurve(V; extent=1, direction=:positive, n=Vec{3, Float64}(0.0,0.0,-3.0), close_loop=false,face_type=:quad)
elseif testCase==2 
    # An imported STL based geometry of an engineering part with various flat faces and some known angles (e.g. 0, 45, and 90 degrees)
    fileName_mesh = joinpath(comododir(),"assets","stl","spur_gear_01.stl")
    M = load(fileName_mesh)
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_,_ = mergevertices(F,V)
end

# E = meshedges(F)
# E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)        
# con_E2F = con_edge_face(F,E_uni)
# con_F2E = con_face_edge(F,E_uni,indReverse)   

A,E,con_E2F = edgefaceangles(F,V; deg=true)


## Visualization
linewidth = 3

Lp = .~isnan.(A) .&& abs.(A).>30
Ap = A[Lp]
Ep = E[Lp]

En,Vn = separate_vertices(Ep,V)
An = simplex2vertexdata(En,Ap,Vn)

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Edge angles")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), color=:white, shading = FastShading, transparency=false)#,strokecolor=:black,strokewidth=0.25)
# normalplot(ax1,F,V)

hp_A = wireframe!(ax1,GeometryBasics.Mesh(Vn,En),linewidth=linewidth, transparency=false, color=An,colormap=:Spectral,colorrange = (-180, 180))

Colorbar(fig[1, 2],hp_A, label = "Angles",ticks =-180:20:180)

fig