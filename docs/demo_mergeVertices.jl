using GLMakie
using Gibbon
using GeometryBasics
using FileIO
using Statistics

"""
dfasdf
"""
# Loading a mesh
fileName_mesh = joinpath(gibbonDir(),"assets","stl","stanford_bunny_low.stl")
M1 = load(fileName_mesh)
F1 = faces(M1)
V1 = coordinates(M1)
F1 = toGeometryBasicsSimplices(F1) 
V1 = toGeometryBasicsPoints(V1) 

F2 = deepcopy(F1)
V2 = deepcopy(V1)


F2,V2,ind1,ind2 = mergeVertices(F2,V2; roundVertices=true)
M2 = GeometryBasics.Mesh(V2,F2)



function face2vertex(F,V,C)
    con_V2F = con_vertex_face(F,V)
    CV = (typeof(C))(undef,length(V))
    for q ∈ eachindex(V)
        CV[q] = mean(C[con_V2F[q]])
    end
    return CV
end

N1,VN=meshnormal(F1,V1)
NV1 = face2vertex(F1,V1,N1)

N2,VN=meshnormal(F2,V2)
NV2 = face2vertex(F2,V2,N2)

## Visualisation

fig = Figure(size=(1200,800))

stepRange = 0:1:50
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 0,linewidth=30)

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

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
poly!(ax1,M1,strokewidth=2,color=:white, shading = FastShading)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
poly!(ax2,M2,strokewidth=2,color=:white, shading = FastShading)


fig

    