using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1 # Sphere 
        nSteps = 50
        # Define sphere using distance function
        xr,yr,zr = ntuple(_->range(-1.0,1.0,nSteps),3)
        A = [norm((x,y,z)) for x in xr, y in yr, z in zr]
        level = 0.5
        cap = false
    elseif testCase == 2 # Torus
        nSteps = 50
        rt = 0.5 # Torus major radius
        # Define torus using distance function 
        n = 100 # Number of points used for distance calculation 
        Vc = circlepoints(rt,n)    
        xr,yr,zr = ntuple(_->range(-1.0,1.0,nSteps),3)    
        A = Array{Float64,3}(undef,(nSteps,nSteps,nSteps))
        for i in 1:1:nSteps
            for j in 1:1:nSteps
                for k in 1:1:nSteps
                    A[i,j,k] = mindist([Point{3,Float64}(xr[i],yr[j],zr[k])],Vc)[1]                
                end
            end
        end
        level = 0.25
        cap = false
    elseif testCase == 3 # Gyroid not capped 
        nSteps = 50
        np = 3 # Number of "periods"
        xr,yr,zr = ntuple(_->range(0,2*pi*np,nSteps),3)
        A = [triplyperiodicminimal((x,y,z), :G) for x in xr, y in yr, z in zr]
        level = 0.0
        cap = false
    elseif testCase == 4 # Gyroid capped    
        nSteps = 50
        np = 3 # Number of "periods"
        xr,yr,zr = ntuple(_->range(0,2*pi*np,nSteps),3)
        A = [triplyperiodicminimal((x,y,z), :G) for x in xr, y in yr, z in zr]
        level = 0.0
        cap = true
    end

    F1,V1 = getisosurface(A; x = xr, y = yr, z = zr, level = level, cap = cap, padValue=1e8)      

    # Visualization
    strokewidth = 0.5
    fig = Figure(size=(800,800))

    stepRange1 = range(minimum(A),maximum(A),30)
    hSlider1 = Slider(fig[2, :], range = stepRange1, startvalue = level,linewidth=30)

    titleString = lift(hSlider1.value) do level
        "level = " * string(level) 
    end

    ax1 = AxisGeom(fig[1, 1]; title=titleString, limits=(minimum(xr),maximum(xr),minimum(yr),maximum(yr),minimum(zr),maximum(zr)))
    hp1 = meshplot!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=strokewidth, strokecolor=:black, stroke_depth_shift=-0.001f0)
    # normalplot(ax1,F1,V1)

    on(hSlider1.value) do level
        F1,V1 = getisosurface(A; x=collect(xr), y=collect(yr), z=collect(zr), level=level, cap=cap, padValue=1e8)    
        F1,V1 = separate_vertices(F1,V1)    
        if !isempty(F1)
            hp1[1] = GeometryBasics.Mesh(V1,F1)       
        end
    end

    slidercontrol(hSlider1,ax1)   

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end