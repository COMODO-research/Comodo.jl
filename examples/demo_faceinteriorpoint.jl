using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using FileIO

using Comodo.Statistics
using Comodo.LinearAlgebra

testCase = 2

# Loading a mesh
if testCase == 1    
    fileName_mesh = joinpath(comododir(),"assets","stl","femur.stl")
elseif testCase == 2
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")    
end

M = load(fileName_mesh)
F = [TriangleFace{Int}(f) for f in faces(M)]
V = [Point{3,Float64}(p) for p in coordinates(M)]
F,V  = mergevertices(F,V)

P_in = [faceinteriorpoint(F,V, i; w=0.5) for i in eachindex(F)]

## Visualization
c = RGBA(0.5, 0.5, 0.5,0.25)

markerSize = 10

fig = Figure(size=(800,800))
ax1 = LScene(fig[1,1]) # ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """rayType = :ray, triSide=1""")
hp1 = poly!(ax1,M,color=c, shading = FastShading, transparency=true,strokecolor=:black, strokewidth=0.25)
hp2 = scatter!(ax1,P_in,color=:red,markersize=markerSize)

# scatter!(ax1,ray_origin,color=:green,markersize=markerSize*2)
# scatter!(ax1,P,color=:red,markersize=markerSize)
# lines!(ax1,P,color=:red)
# lines!(ax1,[ray_origin, ray_origin+3*ray_vector],color=:blue)

stepRange = range(0.0,1.0,20)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 0.5,linewidth=30)

on(hSlider.value) do sliderVal 
    hp2[1] = [faceinteriorpoint(F,V, i; w=sliderVal ) for i in eachindex(F)]    
end

fig
