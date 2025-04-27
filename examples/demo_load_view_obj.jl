using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_triangulated.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
    elseif testCase == 2
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_quadrangulated.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
    elseif testCase == 3
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
    elseif testCase == 4
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh_texture.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
    elseif testCase == 5
        fileName_mesh = joinpath(comododir(),"assets","obj","lego_figure.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","lego_figure.png")    
    end

    M = load(fileName_mesh)
    T = load(fileName_texture)

    ## Visualization
    fig = Figure(size=(800,800))
    # ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Spot the cow")
    ax1 = LScene(fig[1,1])
    hp1 = poly!(ax1,M, color=T, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0)
   
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end