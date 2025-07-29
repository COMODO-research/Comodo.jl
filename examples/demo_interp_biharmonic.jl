using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Random

#=
In this demo biharmonic interpolation is used for 3D scattered data 
interpolation. First a set of random 3D points and "value data" is defined. Next
a grid of points are defined onto which this value data is interpolated. 
=#

Random.seed!(1) # Set seed so demo performs the same each time

GLMakie.closeall()
for testCase = 1:2
    if testCase == 1 # SIngle triangle    
        # Define raw data
        m = 150 # Number of input data points
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,m) # Initialize point vector
        C = Vector{Float64}(undef,m) # Initialize value/color vector
        for q = 1:m # Loop over all points
            V[q] = 2.0*pi*rand(3) # Assign a randon point from 0-2*pi 
            C[q] = sin(2*sqrt(sum(V[q].^2))) # An interesting function with distance from origin
        end

        # Define points for interpolation 
        n = 25
        Vi = gridpoints(range(0.0,2.0*pi,n),range(0.0,2.0*pi,n),range(0,2.0*pi,5)) # Define a grid of points 

        # Interpolate the data onto the grid
        Ci = interp_biharmonic(V,C,Vi)

        # Visualize 
        fig = Figure(size = (800, 800))
        ax = AxisGeom(fig[1, 1])
        hs1 = scatter!(ax, V,markersize=35,color=C,colormap = :Spectral,strokewidth=2)
        hs2 = scatter!(ax, Vi,markersize=15,color=Ci,colormap = :Spectral)
        Colorbar(fig[1, 2],hs1,label = "Value data")
        Legend(fig[1, 3],[hs1,hs2],["Raw","Interpolated"])
        fig

    elseif testCase == 2 # n-slice pizza triangle set
        # Define raw data
        fSin(x,f,a) = 2.0*sin(f*x+a) # Function to change shape of loop
        a = 0.0
        f = 1
        r = 3.0
        n = 50 # Number of input data points
        V = circlepoints(r,n)
        C = Vector{Float64}(undef,n)
        for (i,v) in enumerate(V)
            C[i] = fSin(v[1],f,a)
            V[i] = Point{3,Float64}(v[1], v[2], C[i])
        end

        # Define points for interpolation 
        pointSpacing = pointspacingmean(V; close_loop=true)
        VT = (deepcopy(V),)
        R = ([1],)
        P = (pointSpacing)
        Fi,Vi,_ = regiontrimesh(VT,R,P)

        V_2D = [Point{2,Float64}(v[1], v[2]) for v in V]
        Vi_2D = [Point{2,Float64}(v[1], v[2]) for v in Vi]

        Ci = interp_biharmonic(V_2D,C,Vi_2D)

        for (i,v) in enumerate(Vi)
            Vi[i] = Point{3,Float64}(v[1], v[2], Ci[i]) 
        end

        # Visualize 
        cMap = :viridis
        fig = Figure(size = (800, 800))
        ax1 = AxisGeom(fig[1, 1])
        hs1 = scatter!(ax1, V, markersize=10, color=C, colormap = cMap,strokewidth=2, colorrange=(minimum(C),maximum(C)), depth_shift=-0.01f0)        
        hs2 = meshplot!(ax1, Fi,Vi; color=Ci, colormap = cMap, colorrange=(minimum(C),maximum(C)))
        Colorbar(fig[1, 2],hs1,label = "Value data")
        Legend(fig[1, 3],[hs1,hs2],["Raw","Interpolated"])
        stepRange = range(-pi,pi,101)

        hSlider = Slider(fig[2, :], range = stepRange, startvalue = 0.0, linewidth=30)

        on(hSlider.value) do a     

            for (i,v) in enumerate(V_2D)
                C[i] = fSin(v[1],f,a)        
            end
            
            Ci = interp_biharmonic(V_2D,C,Vi_2D)
            Vi2 = deepcopy(Vi)
            for (i,v) in enumerate(Vi2)
                Vi2[i] = Point{3,Float64}(v[1], v[2], Ci[i]) 
            end 

            V2 = deepcopy(V)
            for (i,v) in enumerate(V2)
                V2[i] = Point{3,Float64}(v[1], v[2], C[i]) 
            end 

            hs1[1] = V2
            hs1.color = C

            hs2[1] = GeometryBasics.Mesh(Vi2,Fi)
            hs2.color = Ci
        end

        slidercontrol(hSlider,ax1)

        fig
    end

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end
