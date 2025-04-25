using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.LinearAlgebra

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

GLMakie.closeall()

for testCase = 1:5
    if testCase == 1
        plateDim1 = [20.0,24.0]
        plateElem1 = [11,16]
        orientation1 = :up
        F1,V1 = quadplate(plateDim1,plateElem1; orientation=orientation1)
        V2 = deepcopy(V1)
        p = eltype(V2)(0.0,0.0,15.0)
        V2 = [v+p for v in V2]

        numSteps = 8
        correspondence = :match

        E = fromtomesh!(F1, V1, V2, numSteps; correspondence = correspondence)
    elseif testCase == 2
        plateDim1 = [20.0,24.0]
        plateElem1 = [11,16]
        orientation1 = :up
        F1,V1 = quadplate(plateDim1,plateElem1; orientation=orientation1)
        V2 = deepcopy(V1)
        p = eltype(V2)(0.0,0.0,15.0)
        V2 = [v+p for v in V2]

        numSteps = 8
        correspondence = :match
        # Non input manipulating function test
        E,V1 = fromtomesh(F1, V1, V2, numSteps; correspondence = correspondence)       
    elseif testCase == 3
        plateDim1 = [20.0,24.0]
        plateElem1 = [11,16]
        orientation1 = :up
        F1,V1 = quadplate(plateDim1,plateElem1; orientation=orientation1)
        ind1 = unique(reduce(vcat,F1))
        p = eltype(V1)(0.0,0.0,15.0)
        V2 = [v+p for v in V1[ind1]]

        numSteps = 8
        correspondence = :faces

        E = fromtomesh!(F1, V1, V2, numSteps; correspondence = correspondence)
    elseif testCase == 4
        r = 1.0
        f = 0.4
        pointSpacing = r/10
        t = f*r
        boxDim = fill(2*f*r,3) # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        F1,V1,C1 = quadbox(boxDim,boxEl)        
        V2 = [v.* (r/norm(v)) for v in V1]
        
        numSteps = 5   
        correspondence = :match

        E = fromtomesh!(F1, V1, V2, numSteps; correspondence = correspondence)
    elseif testCase == 5
        r = 1.0
        f = 0.4
        pointSpacing = r/10
        t = f*r
        boxDim = fill(2*f*r,3) # Dimensionsions for the box in each direction
        F1,V1,C1 = tribox(boxDim,pointSpacing)
        V2 = [v.* (r/norm(v)) for v in V1]

        numSteps = 5   
        correspondence = :match

        E = fromtomesh!(F1, V1, V2, numSteps; correspondence = correspondence)
    end

    # append!(E,fromtomesh!(F, V, V2, numSteps))
    
    Fn = element2faces(E)

    # Visualisation
    cmap = Makie.Categorical(:Spectral) 
    a = 0.25
    fig = Figure(size=(800,800))

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Surface sets to mesh from-to")
    hp1 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=3,shading=FastShading, strokecolor=:black, color=:white, transparency=false, overdraw=false)
    hp2 = scatter!(ax1,V2, color = :red, markersize=15)

    ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Volumetric mesh (cut view)")
    hpSlider = Vector{Any}()
    if isa(Fn,Tuple)
        hp31 = poly!(ax2, GeometryBasics.Mesh(V1,Fn[1]), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
        hp32 = poly!(ax2, GeometryBasics.Mesh(V1,Fn[2]), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
        push!(hpSlider,hp31)
        push!(hpSlider,hp32)
    else
        hp3 = poly!(ax2, GeometryBasics.Mesh(V1,Fn), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)
        push!(hpSlider,hp3)
    end

    pointSpacing = mean(edgelengths(F1,V1))
    VE  = simplexcenter(E,V1)
    ZE = [v[3] for v in VE]
    Z = [v[3] for v in V1]
    zMax = maximum(Z)
    zMin = minimum(Z)
    numSlicerSteps = 3*ceil(Int,(zMax-zMin)/pointSpacing)

    stepRange = range(zMin,zMax,numSlicerSteps)
    hSlider = Slider(fig[2, :], range = stepRange, startvalue = stepRange[end],linewidth=30)

    on(hSlider.value) do z     
        indShow = findall(ZE .<= z)
        if isempty(indShow)
            for hp in hpSlider
                hp.visible=false        
            end
        else
            Fn = element2faces(E[indShow])        
            if isa(Fn,Tuple)
                for hp in hpSlider
                    hp.visible=true    
                end
                Fns,Vns = separate_vertices(Fn[1],V1)
                hpSlider[1][1] = GeometryBasics.Mesh(Vns,Fns)
                Fns,Vns = separate_vertices(Fn[2],V1)
                hpSlider[2][1] = GeometryBasics.Mesh(Vns,Fns)
            else
                hpSlider[1][1].visible=true
                Fn = element2faces(E[indShow])        
                Fns,Vns = separate_vertices(Fn,V1)
                hpSlider[1][1] = GeometryBasics.Mesh(Vns,Fns)
            end
        end
    end

    slidercontrol(hSlider,ax2)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end