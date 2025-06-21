using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.GLMakie.FileIO
using Random

Random.seed!(1) # Set seed so demo performs the same each time

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1 # Triangle mesh bunny
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)
        F = faces(M)
        V = coordinates(M)    
        V = topoints(V) 
        F = tofaces(F)
        F,V,_,_ = mergevertices(F,V)
    elseif testCase == 2 # Quad mesh sphere
        F,V = subquadsphere(4,100)
    end

    V0 = deepcopy(V)

    V = map(x-> x.+ (5.0 .* randn(eltype(V))),V)

    # Set smoothing parameters
    tol = 1e-3
    nMax = 15 # Maximum number of iterations

    λ = 0.5 # Laplacian smoothing parameter

    # HC parameters
    α = 0.1
    β = 0.5

    # Use constrained nodes for the bottom half
    testMode = 2
    if testMode == 2
        z = [v[3] for v in V]
        constrained_points = findall(z.<mean(z))
    else
        constrained_points = nothing
    end

    ## VISUALISATION

    strokeWidth1 = 1

    Vs_HC = smoothmesh_hc(F,V,nMax,α,β; tolDist=tol, constrained_points = constrained_points)
    Vs_LAP = smoothmesh_laplacian(F,V,nMax,λ; constrained_points = constrained_points)
    Ds = [sqrt(sum((Vs_HC[i]-V0[i]).^2)) for i ∈ eachindex(V)]
    Ds_LAP = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
    cLim = maximum(Ds_LAP).*(0.0,1.0)

    fig = Figure(size=(1200,1200))

    ax1 = AxisGeom(fig[1, 1], title = "Original")
    hp1 = meshplot!(ax1,GeometryBasics.Mesh(V0,F),strokewidth=strokeWidth1)

    ax2 = AxisGeom(fig[1, 2], title = "Noisy")
    hp2 = meshplot!(ax2,GeometryBasics.Mesh(V,F),strokewidth=strokeWidth1)

    ax3 = AxisGeom(fig[2, 1])
    hp3 = meshplot!(ax3,GeometryBasics.Mesh(Vs_LAP,F),strokewidth=strokeWidth1,color=Ds_LAP,colormap=Makie.Reverse(:Spectral),colorrange=cLim)

    ax4 = AxisGeom(fig[2, 2])
    hp4 = meshplot!(ax4,GeometryBasics.Mesh(Vs_HC,F),strokewidth=strokeWidth1,color=Ds, colormap=Makie.Reverse(:Spectral),colorrange=cLim)

    Colorbar(fig[:, 3],hp3,label = "Distance")

    stepRange = 0:1:nMax
    hSlider = Slider(fig[3, :], range = stepRange, startvalue = nMax,linewidth=30)

    slidercontrol(hSlider,fig)

    on(hSlider.value) do stepIndex
        # Update first plot
        Vs_LAP = smoothmesh_laplacian(F,V,stepIndex,λ; constrained_points = constrained_points)
        Ds_LAP = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
        ax3.title= "Laplacian smoothed n = " * string(stepIndex) * " times"
        hp3.color = Ds_LAP
        hp3[1] = GeometryBasics.Mesh(Vs_LAP,F)

        # Update second plot
        Vs_HC = smoothmesh_hc(F,V,stepIndex,α,β; tolDist=tol, constrained_points = constrained_points)
        Ds_HC = [sqrt(sum((Vs_HC[i]-V0[i]).^2)) for i ∈ eachindex(V)]
        ax4.title= "HC smoothed n = " * string(stepIndex) * " times"
        hp4.color = Ds_HC
        hp4[1] = GeometryBasics.Mesh(Vs_HC,F)
    end
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end