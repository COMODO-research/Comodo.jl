using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra
using Statistics

# Example geometry
testCase = 1
if testCase == 1 # Triangulated sphere
    r = 1
    F,V = geosphere(3,r)  
    println("Theoretical volume: " * string(4/3*pi*r^3))    
elseif testCase==2 # quad sphere 
    r = 1.0
    F,V = quadsphere(5,r)    
    println("Theoretical volume: " * string(4/3*pi*r^3))    
elseif testCase==3 # quad cube
    r=2*sqrt(3)/2
    F,V = cube(r)    
    println("Theoretical volume: " * string((2*(r./(sqrt(3))))^3))
elseif testCase==4 # triangulated cube
    r=2*sqrt(3)/2
    F,V = cube(r)    
    F = quad2tri(F,V; convert_method = :angle)
    println("Theoretical volume: " * string((2*(r./(sqrt(3))))^3))
elseif testCase==5 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_ = mergevertices(F,V)
end


vol = surfacevolume(F,V)

println("Computed volume: " *string(vol))

# M = GeometryBasics.Mesh(V,F)

# A = facearea(F,V)
# println("Total mesh area: " *string(sum(A)))

# Fn,Vn = separate_vertices(F,V)

# ## Visualization

# fig = Figure(size=(800,800))

# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Face area")
# hp1 = poly!(ax1,GeometryBasics.Mesh(Vn,Fn), color=simplex2vertexdata(Fn,A), shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
# Colorbar(fig[1, 2],hp1)
# fig


