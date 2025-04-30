using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1    
        V = Point{3,Float64}[ [0.0,0.0,0.0], [10.0,0.0,0.0], [10.0,10.0,0.0]]
        close_loop = false
    elseif testCase == 2
        V = Point{3,Float64}[ [0.0,0.0,0.0], [10.0,0.0,0.0], [10.0,10.0,0.0], [0.0,10.0,0.0]]
        close_loop = false
    elseif testCase == 3
        V = Point{3,Float64}[ [0.0,0.0,0.0], [10.0,0.0,0.0], [10.0,10.0,0.0], [0.0,10.0,0.0]]
        close_loop = true
    elseif testCase == 4
        V = Point{3,Float64}[ [0.0,0.0,0.0], [1.0,0.0,0.0], [1.0,2.0,0.0], [2.0,2.0,0.0]]
        close_loop = false
    end

    rMax = 0.5 
    n = 20
    h = 2.0
    VC = filletcurve(V; rMax=rMax,  constrain_method = :max, n=n, close_loop = close_loop, eps_level = 1e-6)
    VC = evenly_sample(VC,50; spline_order = 2)

    Fe,Ve = extrudecurve(VC; extent=h, direction=:both, close_loop=false,face_type=:quad)
    
    # Visualisation
    fig = Figure(size=(1000,1000))
    ax1 = Axis3(fig[1, 1],aspect = :data,title="Filleting a curve")

    hp11 = lines!(ax1, V,linewidth=2,color=:black)
    hp12 = scatter!(ax1, V,markersize=15,color=:black)

    Fes,Ves = separate_vertices(Fe,Ve)
    hp3 = poly!(ax1,GeometryBasics.Mesh(Ves,Fes), strokewidth=1,color=:lightgreen, shading = FastShading,transparency=false)
    # # scatter!(ax1, V[1],markersize=25,color=:yellow)
    # # scatter!(ax1, V[end],markersize=25,color=:red)

    indPlot = collect(1:length(VC))
    if close_loop == true
        push!(indPlot,1)
    end
    hp2 = lines!(ax1, VC[indPlot],linewidth=6,color=:blue)

    stepRange1 = range(0,1.5,500)
    hSlider1 = Slider(fig[2, 1], range = stepRange1, startvalue = 0,linewidth=30)

    on(hSlider1.value) do stepIndex1
        VC = filletcurve(V; rMax=stepIndex1,  constrain_method = :max, n=n, close_loop = close_loop, eps_level = 1e-6)
        VC = evenly_sample(VC,53; spline_order = 2)
        indPlot = collect(1:length(VC))
        if close_loop == true
            push!(indPlot,1)
        end
        hp2[1] = VC[indPlot]

        Fe,Ve = extrudecurve(VC; extent=h, direction=:both, close_loop=false,face_type=:quad)
        Fes,Ves = separate_vertices(Fe,Ve)
        hp3[1] = GeometryBasics.Mesh(Ves,Fes)
    end

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end