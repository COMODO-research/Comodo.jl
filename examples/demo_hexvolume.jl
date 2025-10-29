using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO

#=
This demo shows the use of `hexvol` to compute the volume for hexahedral 
elements. 
=#

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1
        r = 2.5
        pointSpacing = 0.5
        E, V = hexsphere(r,pointSpacing)   
        println("Theoretical volume             : " * string(4/3*pi*r^3))      
    elseif testCase == 2
        r1 = 2.0
        r2 = r1/2
        pointSpacing = 0.25
        E, V = hexspherehollow(r1, r2, pointSpacing)
        println("Theoretical volume             : " * string(4/3*pi*r1^3 - 4/3*pi*r2^3))
    elseif testCase == 3
        pointSpacing = 0.5
        boxDim = [2.5,3.1,4] # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)
        println("Theoretical volume             : " * string(prod(boxDim)))
    elseif testCase == 4
        h = 1.0
        n = 4
        r = 2.5
        nh = 3
        E, V, F, Fb, Cb = hexcylinder(r, h, n; nh=nh, direction=:both)  
        println("Theoretical volume             : " * string(π*r^2*h))
    end

    vol = hexvolume(E,V)

    println("Computed volume from hexahedra      : " *string(sum(vol)))
    println("Mean volume from hexahedra          : " *string(mean(vol)))
    println("---------------------------------------")
    ## Visualization

    cmap = cgrad(:Spectral, 250)

    F = element2faces(E) # Triangular faces
    CE_F = repeat(vol,inner=6)

    Fs,Vs = separate_vertices(F,V)
    CE_Vs = simplex2vertexdata(Fs,CE_F)
    M = GeometryBasics.Mesh(Vs,Fs)

    strokewidth = 1 

    fig = Figure(size=(1200,1200))

    ax1 = AxisGeom(fig[1, 1], title = "Cut mesh")
    hp2 = meshplot!(ax1, Fs, Vs, color=CE_Vs, strokewidth=strokewidth, colorrange = (0,maximum(vol)), colormap=cmap)

    VE  = simplexcenter(E,V)
    ZE = [v[3] for v in VE]
    Z = [v[3] for v in V]
    zMax = maximum(Z)
    zMin = minimum(Z)
    numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F,V)))

    stepRange = range(zMin,zMax,numSlicerSteps)    
    hSlider = Slider(fig[2, 1], range = stepRange, startvalue = numSlicerSteps,linewidth=30)

    Colorbar(fig[1, 2], hp2, ticks = range(0,maximum(vol),25))

    on(hSlider.value) do z 

        B = ZE .<= z
        indShow = findall(B)
        if isempty(indShow)
            hp2.visible=false        
        else        
            hp2.visible=true
            Fs = element2faces(E[indShow])
            Cs = repeat(vol[indShow], inner=6)

            indB = boundaryfaceindices(Fs)        
            Fs = Fs[indB]
            Cs = Cs[indB]
            Fs,Vs = separate_vertices(Fs,V)
            CE_Vs = simplex2vertexdata(Fs,Cs)

            Ms = GeometryBasics.Mesh(Vs,Fs)
            hp2[1] = Ms
            hp2.color = CE_Vs
        end

    end    
    slidercontrol(hSlider,ax1)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end