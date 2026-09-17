using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.LinearAlgebra

## Define example input
F1, V1 = platonicsolid(2, 1.0) # Get an example quadrilateral mesh (for a cube in this case)

fileName_mesh = joinpath(comododir(), "assets", "obj", "spot_quadrangulated.obj")

function import_obj(loadName)
    file_IO = open(loadName)
    F = Vector{NgonFace}()
    V = Vector{Point{3,Float64}}()
    for line in eachline(file_IO)        
        splitLine = strip.(split(line, " "))
        if splitLine[1] == "f"            
            f = Vector{Int}()
            for s in splitLine[2:end] 
                s_split = split(s, "/")
                push!(f, parse(Int64, s_split[1]))
            end
            push!(F, NgonFace{length(f), Int}(f))
        elseif splitLine[1] == "v"
            push!(V, Point{3,Float64}(parse.(Float64, splitLine[2:end])))
        end
    end
    close(file_IO)
    return F, V
end

F2, V2 = import_obj(fileName_mesh)
F2 = [f for f in F2]
V2 = [Point{3, Float64}(p[3], p[1], p[2]) for p in V2]

n = 1

F1n, V1n = subquad(F1, V1, n; method=:Catmull_Clark)
F2n, V2n = subquad(F2, V2, n; method=:Catmull_Clark)

# Visualisation
GLMakie.closeall()

fig = Figure(size=(1200, 800))
ax1 = AxisGeom(fig[1, 1], title="n = " * string(n) * " refinement steps")
hp1 = meshplot!(ax1, F1n, V1n; strokewidth=0.5)
hp2 = edgeplot!(ax1, F1, V1, linewidth=3, color=:red, depth_shift=-0.01f0)


ax2 = AxisGeom(fig[1, 2], title="n = " * string(n) * " refinement steps")
hp3 = meshplot!(ax2, F2n, V2n; strokewidth=0.5)
hp4 = edgeplot!(ax2, F2, V2, linewidth=1, color=:red, depth_shift=-0.01f0)

stepRange = 0:1:3
hSlider = Slider(fig[2, :], range=stepRange, startvalue=1, linewidth=30)

on(hSlider.value) do n
    F1n, V1n = subquad(F1, V1, n; method=:Catmull_Clark)
    F2n, V2n = subquad(F2, V2, n; method=:Catmull_Clark)

    ax1.title = "n = " * string(n) * " refinement steps"
    ax2.title = "n = " * string(n) * " refinement steps"
    hp1[1] = GeometryBasics.Mesh(V1n, F1n)
    hp3[1] = GeometryBasics.Mesh(V2n, F2n)
end

fig