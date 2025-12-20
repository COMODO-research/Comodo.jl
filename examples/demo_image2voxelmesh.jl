using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1
        nSteps = 25
        w = 2.0
        xr,yr,zr = ntuple(_->range(-w, w,nSteps),3)
        I = [sqrt(x^2 + y^2 + z^2) for x in xr, y in yr, z in zr]    
        d = 2.0
        voxelSelection = [i<=d for i in I] # Bool array
    elseif testCase == 2
        nSteps = 75
        w = 4.77
        xr,yr,zr = ntuple(_->range(-w, w,nSteps),3)
        function thisImage(x,y,z)
            phi=(1+sqrt(5))/2
            return 1.0/6.0 *(2.0 - (cos(x + phi*y) + cos(x - phi*y) + cos(y + phi*z) + cos(y - phi*z) + cos(z - phi*x) + cos(z + phi*x)))        
        end
        I = [thisImage(x,y,z) for x in xr, y in yr, z in zr]    
        B = [i<0.0 for i in I]
        voxelSelection = findall(B) # Cartesian indices
   end

    meshType = :boundaryfaces
    voxelSize = ((2.0*w)/(nSteps-1), (2.0*w)/(nSteps-1), (2.0*w)/(nSteps-1))
    origin = Point{3,Float64}(-w, -w, -w)

    M, V, CM = image2voxelmesh(I, voxelSelection; meshType=meshType, voxelSize=voxelSize, origin=origin)

    # Visualization
    if meshType == :elements
        F = element2faces(M)
        CF = repeat(CM, inner=6)
    else # meshType == :faces || meshType == :boundaryfaces
        F = M
        CF = CM
    end
    Fs, Vs = separate_vertices(F,V)
    Cs = simplex2vertexdata(Fs,CF)

    fig = Figure(size=(800,800))   
    ax1 = AxisGeom(fig[1,1]; limits=(-w-voxelSize[2],w+voxelSize[2],-w-voxelSize[1],w+voxelSize[1],-w-voxelSize[3],w+voxelSize[3]))
    hp1 = meshplot!(ax1, Fs, Vs, color=Cs, strokewidth=1.0, colormap=Reverse(:Spectral))
    # scatter!(ax1,V, markersize=10, color=:red)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end