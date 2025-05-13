using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using FileIO

c1 = RGBf(1.0, 0.30196078431372547, 0.023529411764705882)
c2 = RGBf(0.2235294117647059, 1.0, 0.0784313725490196)
c3 = RGBA(0.5803921568627451, 0.3411764705882353, 0.9215686274509803,0.9)

#=
This demo shows the use of the dualclag function. 
=#

GLMakie.closeall()

for testCase = 1:7
    if testCase == 1 
        F,V = geosphere(2,1.0)
        con_type = :face
    elseif testCase == 2 
        F,V = geosphere(2,1.0)
        con_type = :edge
    elseif testCase == 3 
        F,V = geosphere(1,1.0)
        B = [v[3]>0 for v in V]
        BF = [all(B[f]) for f in F]
        F = F[BF]
        F,V = remove_unused_vertices(F,V)
        con_type = :face
    elseif testCase == 4 
        F,V = geosphere(1,1.0)
        B = [v[3]>0 for v in V]
        BF = [all(B[f]) for f in F]
        F = F[BF]
        F,V = remove_unused_vertices(F,V)
        con_type = :edge
    elseif testCase == 5
        F,V = tetrahedron(√3)     
        con_type = :face   
    elseif testCase == 6
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V = mergevertices(F,V)
        # n = 1
        # F,V = subtri(F,V,n; method=:Loop)
        con_type = :face
    elseif testCase == 7
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V = mergevertices(F,V)
        B = [v[3]>-10 for v in V]
        BF = [all(B[f]) for f in F]
        F = F[BF]
        F,V = remove_unused_vertices(F,V)    
        con_type = :face
    # elseif testCase == 8
    #     # Loading a mesh
    #     fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    #     M = load(fileName_mesh)

    #     # Obtain mesh faces and vertices
    #     F = tofaces(faces(M))
    #     V = topoints(coordinates(M))
    #     F,V,_ = mergevertices(F,V)
    #     con_type = :face
    end

    s = 0.5
    Fs,Fq,Vs = dualclad(F,V,s; connectivity=con_type)

    # Visualisation
    strokewidth = 0.5

    fig = Figure(size = (1200,1200))
    ax1 = AxisGeom(fig[1, 1], title = "dual clad surface")
    hp1 = meshplot!(ax1, F, V; color=c3, strokewidth=0, transparency=true)
    hp2 = meshplot!(ax1, Fs, Vs; color=c1, strokewidth=strokewidth, depth_shift=-0.01f0)
    hp3 = meshplot!(ax1, Fq, Vs; color=c2, strokewidth=strokewidth, depth_shift=-0.01f0)

    stepRange = range(1.0,0.0,50)
    hSlider = Slider(fig[2, :], range = stepRange, startvalue = s,linewidth=30)

    slidercontrol(hSlider,fig)

    on(hSlider.value) do s
        Fs,Fq,Vs = dualclad(F,V,s; connectivity=con_type)        
        hp2[1] = GeometryBasics.Mesh(Vs,Fs)
        hp3[1] = GeometryBasics.Mesh(Vs,Fq)
    end

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end


