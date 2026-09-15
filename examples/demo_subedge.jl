using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.GLMakie.Colors

using FileIO


GLMakie.closeall()

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        np = 10
        E = [LineFace{Int}(i, i+1) for i in 1:np-1]
        V = [Point{3, Float64}(cos(t)-1.0, sin(t), t/(2*pi)-0.5) for t in range(0.0, 2*pi, np)]
    elseif testCase == 2
        np = 10
        E = [LineFace{Int}(i, i+1) for i in 1:np-1]
        push!(E, LineFace{Int}(1,np))
        V = [Point{3, Float64}(cos(t)-1.0, sin(t), t/(2*pi)-0.5) for t in range(0.0, 2*pi, np)]
    elseif testCase==3
        V = [   Point{3, Float64}( 0.0,  0.0,  0.0),
                Point{3, Float64}( 0.5,  0.0,  0.0), 
                Point{3, Float64}( 1.0,  0.0, -0.5), 
                Point{3, Float64}( 0.0,  0.5,  0.0),
                Point{3, Float64}( 0.0,  1.0,  0.5),
                Point{3, Float64}(-0.5,  0.0,  0.0),
                Point{3, Float64}(-1.0,  0.0,  0.5),
                Point{3, Float64}( 0.0, -0.5,  0.0),
                Point{3, Float64}( 0.0, -1.0, -0.5),
                ]

        E = [   LineFace{Int}(1, 2),
                LineFace{Int}(2, 3),
                LineFace{Int}(1, 4),
                LineFace{Int}(4, 5),
                LineFace{Int}(1, 6),
                LineFace{Int}(6, 7), 
                LineFace{Int}(1, 8),
                LineFace{Int}(8, 9), 
                ]
    end

    n = 2

    Es1, Vs1 = subedge(E, V, n; method=:linear)
    Es2, Vs2 = subedge(E, V, n; method=:smooth)

    fig1 = Figure(size = (1200,800))

    ax1 = AxisGeom(fig1[1, 1], title="Linear splitting, n=$n")
    scatter!(ax1, V, color=:blue, markersize=20)
    h11 = scatter!(ax1, Vs1, color=:red, markersize=12)
    h12 = edgeplot!(ax1, Es1, Vs1, color=:black, linewidth=3)

    ax2 = AxisGeom(fig1[1, 2], title="Smooth splitting, n=$n")
    scatter!(ax2, V, color=:blue, markersize=20)
    h21 = scatter!(ax2, Vs2, color=:red, markersize=12)
    h22 = edgeplot!(ax2, Es2, Vs2, color=:black, linewidth=3)


    stepRange1 = 0:6
    hSlider1 = Slider(fig1[2, :], range = stepRange1, startvalue = n, linewidth=30)

    on(hSlider1.value) do n
        Es1, Vs1 = subedge(E, V, n; method=:linear)
        Es2, Vs2 = subedge(E, V, n; method=:smooth)

        h11[1] = Vs1
        h12[1] = GeometryBasics.Mesh(Vs1, Es1)
        ax1.title = "Linear splitting, n=$n"

        h21[1] = Vs2
        h22[1] = GeometryBasics.Mesh(Vs2, Es2)
        ax2.title = "Smooth splitting, n=$n"
    end

    screen = display(GLMakie.Screen(), fig1)

end
