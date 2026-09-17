using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

F1, V1 = platonicsolid(4, 1.0)

fileName_mesh = joinpath(comododir(), "assets", "obj", "spot_control_mesh.obj")
M = load(fileName_mesh)
# Obtain mesh faces and vertices
F2 = tofaces(faces(M))
V2 = [Point{3, Float64}(p[3], p[1], p[2]) for p in coordinates(M)]
F2, V2, _, _ = mergevertices(F2, V2)

n = 1

F1n, V1n = subtri(F1, V1, n; method=:Loop)
F2n, V2n = subtri(F2, V2, n; method=:Loop)

# Visualisation
GLMakie.closeall()

fig = Figure(size=(1200, 800))
ax1 = AxisGeom(fig[1, 1], title="n = " * string(n) * " refinement steps")
hp1 = meshplot!(ax1, F1n, V1n; strokewidth=0.5)
hp2 = edgeplot!(ax1, F1, V1, linewidth=3, color=:red, depth_shift=-0.01f0)


ax2 = AxisGeom(fig[1, 2], title="n = " * string(n) * " refinement steps")
hp3 = meshplot!(ax2, F2n, V2n; strokewidth=0.25)
hp4 = edgeplot!(ax2, F2, V2, linewidth=1, color=:red, depth_shift=-0.01f0)

stepRange = 0:1:3
hSlider = Slider(fig[2, :], range=stepRange, startvalue=1, linewidth=30)

on(hSlider.value) do n
    F1n, V1n = subtri(F1, V1, n; method=:Loop)
    F2n, V2n = subtri(F2, V2, n; method=:Loop)
    
    ax1.title = "n = " * string(n) * " refinement steps"
    ax2.title = "n = " * string(n) * " refinement steps"
    hp1[1] = GeometryBasics.Mesh(V1n, F1n)
    hp3[1] = GeometryBasics.Mesh(V2n, F2n)
end

fig