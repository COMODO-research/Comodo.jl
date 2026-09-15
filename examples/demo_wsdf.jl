using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Rotations

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1 # Open cylinder 
        # Upward line point set
        n = 20
        h = 15.0
        V1 = collect(range(Point{3,Float64}(0.0, 0.0, -h/2.0), Point{3,Float64}(0.0, 0.0, h/2.0), n))
        r = 10.0
        R1 = fill(10.0, n)

        # Define grid ranges
        voxelSize = (1.0, 1.0, 1.0)
        xr = -1.2*r:voxelSize[1]:1.2*r
        yr = -1.3*r:voxelSize[2]:1.3*r
        zr = -h/2.0:voxelSize[3]:h/2.0

        rMax = maximum(R1)
    elseif testCase == 2 # Doughnut
        voxelSize = (1.0, 1.0, 1.0)

        # Create example curve 
        r1 = 30.0 # Radius
        n1 = ceil(Int, (2π*r1)/minimum(voxelSize))+1 # Number of points
        V1 = circlepoints(r1,n1; dir=:acw)

        # Create example normalisation factors or radii
        R1 = fill(20.0, length(V1))

        # Define grid for distance computation, here based on input curve
        rMax = maximum(R1)
        p_min = minp(V1)
        p_max = maxp(V1)
        numVoxelsAdd = 2
        p_offset = Point{3,Float64}(rMax+numVoxelsAdd*voxelSize[1], rMax+numVoxelsAdd*voxelSize[2], rMax+numVoxelsAdd*voxelSize[3])
        p_origin = p_min - p_offset
        p_end = p_max + p_offset

        # Define grid ranges
        xr = p_origin[1]:voxelSize[1]:p_end[1]
        yr = p_origin[2]:voxelSize[2]:p_end[2]
        zr = p_origin[3]:voxelSize[3]:p_end[3]
    elseif testCase == 3 # Funcky doughnut
        voxelSize = (1.0, 1.0, 1.0)

        # Create example curve 
        r1 = 30.0 # Radius
        n1 = ceil(Int, (2π*r1)/minimum(voxelSize))+1 # Number of points
        V1 = circlepoints(r1,n1; dir=:acw)
        rMin = 5.0
        rMax = 15.0
        t = range(0.0, 2π, n1)

        # Create example normalisation factors or radii
        rFluc = rMax-rMin
        f = 3.0 # frequency
        R1 = rMin .+ rFluc/2.0 .+ rFluc/2.0.*cos.(f*t)

        # Define grid for distance computation, here based on input curve
        rMax = maximum(R1)
        p_min = minp(V1)
        p_max = maxp(V1)
        numVoxelsAdd = 2
        p_offset = Point{3,Float64}(rMax+numVoxelsAdd*voxelSize[1], rMax+numVoxelsAdd*voxelSize[2], rMax+numVoxelsAdd*voxelSize[3])
        p_origin = p_min - p_offset
        p_end = p_max + p_offset

        # Define grid ranges
        xr = p_origin[1]:voxelSize[1]:p_end[1]
        yr = p_origin[2]:voxelSize[2]:p_end[2]
        zr = p_origin[3]:voxelSize[3]:p_end[3]
    elseif testCase == 4 # Aorta like
        voxelSize = (2.0, 2.0, 2.0)

        pointSpacing = 2.0 
        rb = 60.0
        L = 200.0
        nc = 20
        V11 = [Point{3, Float64}(rb*cos(t), rb*sin(t), 0.0) for t in range(pi, 0.0, nc)]
        push!(V11, Point{3, Float64}(rb, -L, 0.0))
        V11 = evenly_space(V11, pointSpacing; close_loop=false, must_points=[1, nc, nc+1])
        R11 = collect(range(30.0, 20.0, length(V11)))

        a2 = ((105.0)/180)*pi
        g2 = ((110.0)/180)*pi
        L2 = 95.0
        V22_1 = Point{3, Float64}(rb*cos(a2), rb*sin(a2), 0.0)
        _, indMin = mindist([V22_1], V11; getIndex=Val(true))
        V22_1 = V11[indMin[1]]
        V22 = [V22_1, V22_1 + Point{3, Float64}(L2*cos(g2), L2*sin(g2), 0.0)]
        V22 = evenly_space(V22, pointSpacing; close_loop=false, spline_order=2)
        R22 = collect(range(10.0, 8.0, length(V22)))
        R22[1] = R11[indMin[1]]

        a3 = ((95.0)/180)*pi
        g3 = ((90.0)/180)*pi
        L3 = 95.0
        V33_1 = Point{3, Float64}(rb*cos(a3), rb*sin(a3), 0.0)
        _, indMin = mindist([V33_1], V11; getIndex=Val(true))
        V33_1 = V11[indMin[1]]
        V33 = [V33_1, V33_1 + Point{3, Float64}(L3*cos(g3), L3*sin(g3), 0.0)]
        V33 = evenly_space(V33, pointSpacing; close_loop=false, spline_order=2)
        R33 = collect(range(10.0, 8.0, length(V33)))
        R33[1] = R11[indMin[1]]

        a4 = ((85.0)/180)*pi
        g4 = ((65.0)/180)*pi
        L4 = 95.0
        V44_1 = Point{3, Float64}(rb*cos(a4), rb*sin(a4), 0.0)
        _, indMin = mindist([V44_1], V11; getIndex=Val(true))
        V44_1 = V11[indMin[1]]
        V44 = [V44_1, V44_1 + Point{3, Float64}(L4*cos(g4), L4*sin(g4), 0.0)]
        V44 = evenly_space(V44, pointSpacing; close_loop=false, spline_order=2)
        R44 = collect(range(10.0, 8.0, length(V44)))
        R44[1] = R11[indMin[1]]

        V1 = deepcopy(V11)
        append!(V1, V22)
        append!(V1, V33)
        append!(V1, V44)

        R1 = deepcopy(R11)
        append!(R1, R22)
        append!(R1, R33)
        append!(R1, R44)

        # Define grid for distance computation, here based on input curve
        rMax = maximum(R1)
        p_min = minp(V1)
        p_max = maxp(V1)
        numVoxelsAdd = 2
        p_offset = Point{3,Float64}(rMax+numVoxelsAdd*voxelSize[1], rMax+numVoxelsAdd*voxelSize[2], rMax+numVoxelsAdd*voxelSize[3])
        p_origin = p_min - p_offset
        p_end = p_max + p_offset

        # Define grid ranges
        xr = p_origin[1]:voxelSize[1]:p_end[1]
        yr = p_origin[2]:voxelSize[2]:p_end[2]
        zr = p_origin[3]:voxelSize[3]:p_end[3]
    end

    # Now computed radius weighted signed distance field
    M = wsdf(xr, yr, zr, V1, R1; closest_type=:weighted)
    siz = size(M)

    # Compute level set surface for visualisation
    FM, VM = getisosurface(M; x=xr, y=yr, z=zr, level=0.0, cap=false, padValue=1e8)      
    # VM = smoothmesh_hc(FM, VM, 25)

    ## Visualization
    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1], title = "Graph with radial data")
    hl = scatter!(ax1, V1, markersize = 10, color = R1, colorrange=(0.0, rMax), colormap=:viridis)
    hm1 = meshplot!(ax1, FM, VM, color=(:white, 0.5), strokewidth=0.0, transparency= true)
    Colorbar(fig[1, 2], hl)

    ax2 = AxisGeom(fig[1, 3], title = "Levelset image")
    sliceIndexMid = round(Int,siz[3]/2.0)
    sliceIndex = sliceIndexMid
    hM = heatmap!(ax2, xr, yr, M[:,:, sliceIndex], interpolate=false, colorrange=(-1.0, 1.0), colormap=:redsblues)
    Makie.translate!(hM, 0.0, 0.0, (sliceIndex-sliceIndexMid)*voxelSize[3])
    hm2 = meshplot!(ax2, FM, VM, color=(:white, 0.5), strokewidth=0.0, transparency= true)
    Colorbar(fig[1, 4], hM)

    hSlider1 = Slider(fig[2, :], range = 1:1:siz[3], startvalue = sliceIndex, linewidth=30)

    on(hSlider1.value) do sliceIndex
        hM[3] = M[:,:,sliceIndex]
        Makie.translate!(hM, 0.0, 0.0, (sliceIndex-sliceIndexMid)*voxelSize[3])
    end

    screen = display(GLMakie.Screen(), fig)

end