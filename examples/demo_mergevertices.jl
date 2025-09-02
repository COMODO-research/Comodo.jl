using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

# Demonstration of the mergevertices function

# Loading a mesh
fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M1 = load(fileName_mesh)
F1 = tofaces(faces(M1)) 
V1 = topoints(coordinates(M1)) 

F2 = deepcopy(F1)
V2 = deepcopy(V1)


V2,indUnique,indMap = mergevertices(V2; pointSpacing=pointspacingmean(F2,V2))
indexmap!(F2,indMap) #Update indices in F2

M2 = GeometryBasics.Mesh(V2,F2)

NV1 = vertexnormal(F1,V1)
NV2 = vertexnormal(F2,V2)

## Visualisation
GLMakie.closeall()

fig = Figure(size=(1200,800))

stepRange = 0:1:50
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 0,linewidth=30)

titleString = lift(hSlider.value) do stepIndex
    "V = V + N*" * string(stepIndex) 
end

M1 = lift(hSlider.value) do stepIndex
    V1n = [V1[i].+ stepIndex.*NV1[i] for i ∈ eachindex(V1)]
    return GeometryBasics.Mesh(V1n,F1)
end

M2 = lift(hSlider.value) do stepIndex
    V2n = [V2[i].+ stepIndex.*NV2[i] for i ∈ eachindex(V2)]
    return GeometryBasics.Mesh(V2n,F2)
end

ax1 = AxisGeom(fig[1, 1],  title = titleString)
poly!(ax1,M1,strokewidth=2,color=:white, shading = FastShading)

ax2 = AxisGeom(fig[1, 2],  title = titleString)
poly!(ax2,M2,strokewidth=2,color=:white, shading = FastShading)

fig[2, :]=hSlider

fig