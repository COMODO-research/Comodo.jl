using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors

#=
This demo shows the use of `tet2hex` to convert tetrahedral elements to hexahedral elements
domain. 
=#

GLMakie.closeall()

for testCase = 1:1
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
    end
    F = element2faces(E)

    Eh,Vh = tet2hex(E,V)
    Fh = element2faces(Eh)

    # Visualisation
    markersize = 25

    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original and converte hexahedral mesh")
    hp1 = meshplot!(ax1, F, V, strokewidth=3, strokecolor=:blue, color=RGBA{Float64}(0.0, 0.0, 1.0, 0.5), transparency = true)
    hp2 = normalplot(ax1, F, V)

    hp1 = meshplot!(ax1, Fh, Vh, strokewidth=2, strokecolor=:red, color=RGBA{Float64}(1.0, 0.0, 0.0, 0.5), transparency = true)
    hp2 = normalplot(ax1, Fh, Vh; color=:red)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end