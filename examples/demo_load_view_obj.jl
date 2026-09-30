using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1 # Mixed mesh
        # fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh_texture.obj")
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = NgonFace                   
        M = load(fileName_mesh; facetype=facetype)
        # M_meta = load(fileName_mesh; facetype=facetype)
        # M = GeometryBasics.Mesh(M_meta.vertex_attributes.position, M_meta.faces; uv=M_meta.vertex_attributes.uv)  
    elseif testCase == 2 # Quads
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_quadrangulated.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = QuadFace{Int}
        M = load(fileName_mesh; facetype=facetype)
    elseif testCase == 3 # Triangles
        fileName_mesh = joinpath(comododir(),"assets","obj","spot_triangulated.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","spot_texture.png")
        facetype = TriangleFace{Int}
        M = load(fileName_mesh; facetype=facetype)
    elseif testCase == 4
        fileName_mesh = joinpath(comododir(),"assets","obj","lego_figure.obj")
        fileName_texture = joinpath(comododir(),"assets","obj","lego_figure.png")    
        facetype = TriangleFace{Int}
        M_meta = load(fileName_mesh; facetype=facetype)
        M = GeometryBasics.Mesh(M_meta.vertex_attributes.position, M_meta.faces; uv=M_meta.vertex_attributes.uv)
    end
    
    T = load(fileName_texture)

    ## Visualization
    fig = Figure(size=(800,800))
    ax1 = AxisGeom(fig[1, 1], title = "Spot the cow")
    
    hp1 = meshplot!(ax1, M, color=T, strokewidth=0.25)
   
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end