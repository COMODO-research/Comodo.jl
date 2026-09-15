using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.Rotations
using FileIO

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        F,V = geosphere(2,1.0)
    elseif testCase == 2
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))        
        F,V,ind1,ind2 = mergevertices(F,V)
        # F,V=subtri(F,V,1)
    elseif testCase == 3
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))    
        F,V,ind1,ind2 = mergevertices(F,V)
    end

    # Input parameters
    p = surface_centroid(F, V); # Point on cutting plane
    n = normalizevector(Vec{3, Float64}(0.0, 1.0, 1.0))# Cutting plane normal
    snapTolerance = 1e-6

    cutType = :full
    Fn, Vn, Cn, En = trisurfslice(F,V,n,p; output_type=cutType)

    ## Visualization
    cmap = colormap = cgrad(:Spectral, 5, categorical = true)
    
    Fns, Vns = separate_vertices(Fn, Vn)
    Cns_V = simplex2vertexdata(Fns, Cn)
    
    s = 1.25*maximum([maximum(map(v-> v[i],V)) - minimum(map(v-> v[i],V)) for i ∈ 1:3])

    R = rotation_between(n,[0.0, 0.0, 1.0])
    plateDim = (s,s)
    plateElem = (1,1)
    FG1,VG1 = quadplate(plateDim,plateElem)
    VGn = [Point{3, Float64}(R'*v)+p for v ∈ VG1]
    
    fig = Figure(size=(1200, 800))

    ax1 = AxisGeom(fig[1, 1], title = "")    
    hp1 = meshplot!(ax1, FG1, VGn; strokewidth=2, strokecolor=:red, color=(:red, 0.1), transparency=true) # Show plate 
    hp2 = meshplot!(ax1, Fns, Vns; color=Cns_V, colorrange = (-2.5, 2.5), colormap=cmap) # Show mesh colored to sliced state
    Colorbar(fig[1, 2], hp2, ticks=-2:1:2)

    ax2 = AxisGeom(fig[1, 3], title = "A sliced mesh")    
    hp3 = meshplot!(ax2, Fn[Cn.<=0], Vn) # Show one side of sliced mesh 
    hp4 = edgeplot!(ax2, En, Vn, linewidth=3, color=:blue) # Show cut edge
    hp5 = meshplot!(ax2, FG1, VGn; strokewidth=2, strokecolor=:red, color=(:red, 0.1), transparency=true) # Show plate 

    Legend(fig[1, 4], [hp3, hp4], ["Sliced surface", "Cut boundary edges"])
    stepRange = range(-s, s, 500)
    hSlider = Slider(fig[2, :], range = stepRange, startvalue = 0,linewidth=30)

    on(hSlider.value) do stepIndex 
        pp = p + stepIndex*n
        Fn, Vn, Cn, En = trisurfslice(F,V,n,pp; output_type=cutType) 
        
        if isempty(Fn)
            Mn = GeometryBasics.Mesh(V,F)
            CnV = zeros(length(V))
        else
            Fns, Vns = separate_vertices(Fn, Vn)
            Cns_V = simplex2vertexdata(Fns, Cn)
            Mn = GeometryBasics.Mesh(Vns, Fns)
        end
        
        VGn = [Point{3, Float64}(R'*v)+pp for v ∈ VG1] # Rotate plane    
        
        hp1[1] = GeometryBasics.Mesh(VGn, FG1)
        hp2[1] = Mn
        hp2.color = Cns_V
        hp3[1] = GeometryBasics.Mesh(Vn, Fn[Cn.<=0])
        if isempty(En)
            hp4.visible=false
        else
            hp4.visible=true
            hp4[1] = GeometryBasics.Mesh(Vn, En)
        end
        hp5[1] = GeometryBasics.Mesh(VGn, FG1)
    end

    slidercontrol(hSlider,ax1)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end