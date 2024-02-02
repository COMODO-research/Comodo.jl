
# Loading a mesh
fileName_mesh = joinpath(gibbonDir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)
F = faces(M)
V = coordinates(M)
F = toGeometryBasicsSimplices(F) 
V = toGeometryBasicsPoints(V) 
F,V,ind1,ind2 = mergeVertices(F,V; roundVertices=true)

V0 = deepcopy(V)

M0 = GeometryBasics.Mesh(V0,F)
V = map(x-> x.+ (10.0 .* rand(eltype(V))),V)
M = GeometryBasics.Mesh(V,F)

E = meshEdges(F)
E_uni,_,_ = unique_simplices(E)
con_V2V = con_vertex_vertex(E_uni)

n = 10
λ = 0.5

function smoothMesh_Laplacian(F,V,λ=0.5, n=1, con_V2V=missing)
    if ismissing(con_V2V)
        E = meshEdges(F)
        E_uni,_,_ = unique_simplices(E)
        con_V2V = con_vertex_vertex(E_uni)
    end
    for _ = 1:1:n
        Vs = deepcopy(V)
        for q ∈ eachindex(V)
            Vs[q] = (1.0-λ).*Vs[q] .+ λ*mean(V[con_V2V[q]])
        end
        V = Vs
    end
    return V
end

function smoothMesh_HC(F,V,α=0.1, β=0.5, n=1, con_V2V=missing)
    if ismissing(con_V2V)
        E = meshEdges(F)
        E_uni,_,_ = unique_simplices(E)
        con_V2V = con_vertex_vertex(E_uni)
    end
    P = deepcopy(V)
    B = deepcopy(V)
    for _ = 1:1:n        
        Q = deepcopy(P)
        for i ∈ eachindex(V)
            P[i] = mean(Q[con_V2V[i]])
            B[i] = P[i] .- (α.*V[i] .+ (1-α).*Q[i])
        end
        
        for i ∈ eachindex(V)            
            P[i] = P[i] .- (β.*B[i] .+ (1-β).* mean(B[con_V2V[i]]))
        end        
    end
    return P
end


α = 0.1
β = 0.5

nMax = 100
# Vs = smoothMesh_Laplacian(F,V,λ,nMax,con_V2V)    
Vs_HC = smoothMesh_HC(F,V,α, β, n, con_V2V)
Vs_LAP = smoothMesh_Laplacian(F,V,λ, n, con_V2V)
Ds = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
cLim = maximum(Ds).*(0.0,1.0)

##############
strokeWidth1 = 1

fig = Figure(size=(900,900))

stepRange = 0:1:nMax
hSlider = Slider(fig[3, :], range = stepRange, startvalue = 0,linewidth=30)

Ms = lift(hSlider.value) do stepIndex
    Vs = smoothMesh_HC(F,V,α, β, stepIndex, con_V2V)
    return GeometryBasics.Mesh(Vs,F)
end

Ms_lap = lift(hSlider.value) do stepIndex
    Vs = smoothMesh_Laplacian(F,V,λ, stepIndex, con_V2V)
    return GeometryBasics.Mesh(Vs,F)
end

Ds = lift(hSlider.value) do stepIndex
    Vs = smoothMesh_HC(F,V,α, β, stepIndex, con_V2V)    
    return [sqrt(sum((Vs[i]-V[i]).^2)) for i ∈ eachindex(V)]
end

Ds_lap = lift(hSlider.value) do stepIndex
    Vs = smoothMesh_Laplacian(F,V,λ, stepIndex, con_V2V)   
    return [sqrt(sum((Vs[i]-V[i]).^2)) for i ∈ eachindex(V)]
end

titleString_lap = lift(hSlider.value) do stepIndex
    "Laplacian smoothed n = " * string(stepIndex) * " times"
end

titleString = lift(hSlider.value) do stepIndex
    "HC smoothed n = " * string(stepIndex) * " times"
end

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original")
poly!(ax1,M0,strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Noisy")
poly!(ax2,M,strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z",title=titleString_lap)
hc1 = poly!(ax3,Ms_lap,strokewidth=strokeWidth1,color=Ds_lap, shading = FastShading,colormap=Makie.Reverse(:Spectral),colorrange=cLim)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z",title=titleString)
hc2 = poly!(ax4,Ms,strokewidth=strokeWidth1,color=Ds, shading = FastShading,colormap=Makie.Reverse(:Spectral),colorrange=cLim)

Colorbar(fig[:, 3],hc1,label = "Distance")

fig[3,:] = hSlider
sliderControl(hSlider,fig)

fig