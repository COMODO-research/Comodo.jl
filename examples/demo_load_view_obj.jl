using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = TriangleFace{Int}                   
    elseif testCase == 2
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_quadrangulated.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = QuadFace{Int}
    elseif testCase == 3
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_triangulated.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = TriangleFace{Int}
    elseif testCase == 4
        fileName_mesh = joinpath(comododir(),"assets","obj","lego_figure.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","lego_figure.png")    
        facetype = TriangleFace{Int}
    elseif testCase == 5 # Mixed mesh 
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh_texture.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = TriangleFace{Int}
    end

    M = load(fileName_mesh; facetype=facetype)
    T = load(fileName_texture)

    ## Visualization
    fig = Figure(size=(800,800))
    ax1 = AxisGeom(fig[1, 1], title = "Spot the cow")
    
    hp1 = poly!(ax1,M, color=T, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1,stroke_depth_shift=-0.001f0)
   
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end