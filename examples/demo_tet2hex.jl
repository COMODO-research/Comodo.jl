using Comodo
using GLMakie
using GeometryBasics

#=
This demo shows the use of `tet2hex` to convert tetrahedral elements to hexahedral elements
domain. 
=#

testCase = 1

if testCase == 1
    E = [Tet4{Int}(1,2,3,4),Tet4{Int}(2,3,4,5),Tet4{Int}(6,7,8,9)]
    V = [Point{3,Float64}(-1.0,0.0,0.0),
         Point{3,Float64}( 1.0,0.0,0.0),
         Point{3,Float64}( 0.0,1.0,0.0),
         Point{3,Float64}( 0.0,0.5,1.0),
         Point{3,Float64}( 1.0,1.0,1.0),
         Point{3,Float64}( 2.0,0.0,0.0),
         Point{3,Float64}( 4.0,0.0,0.0),
         Point{3,Float64}( 3.0,1.0,0.0),
         Point{3,Float64}( 3.0,0.5,1.0),
         ]
elseif testCase == 2
    E = [hex8{Int}(1,2,3,4,5,6,7,8)]
    V = [Point{3,Float64}(0.0,0.0,0.0),
         Point{3,Float64}(1.0,0.0,0.0),
         Point{3,Float64}(1.0,1.0,0.0),
         Point{3,Float64}(0.0,1.0,0.0),         
         Point{3,Float64}(0.0,0.0,1.0),
         Point{3,Float64}(1.0,0.0,1.0),
         Point{3,Float64}(1.0,1.0,1.0),
         Point{3,Float64}(0.0,1.0,1.0),
         ]
end
F = element2faces(E)

Eh,Vh = tet2hex(E,V)
Fh = element2faces(Eh)

# Visualisation
markersize = 25

fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Wireframe of hexahedral mesh")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=true, overdraw=false)
hp2 = normalplot(ax1,F,V)

hp1 = poly!(ax1,GeometryBasics.Mesh(Vh,Fh), strokewidth=2,shading=FastShading,strokecolor=:red, color=:red, transparency=true, overdraw=false)
hp2 = normalplot(ax1,Fh,Vh;color=:red)

fig
