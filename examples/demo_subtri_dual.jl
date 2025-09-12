using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics
using FileIO

#=
This demo shows the use of `subtri_dual` to refine triangulated meshes. The
dual subtriangulation method is sometimes referred to as √3-subdevision (see 
also [1]). Following 2k subdevision steps each original triangle now represents 
9ᵏ triangles. The implementation here largely follows reference [1] however 
smoothing can be turned on/off using the `smooth` option (default is `true`). 
In addition, there is the option to contrain the boundary during smoothing by 
setting `constrain_boundary=true` (default `false`). Finally there is the option 
to split the boundary edges by setting `split_boundary=true` (default is 
`false`).      
# References    
[1] Kobbelt, √3-subdivision, 2000, SIGGRAPH 2000. 
=#

GLMakie.closeall()
for testCase = 1:4
    if testCase == 1
        r = 1.0 # Radius
        n1 = 0
        F,V = tridisc(r,n1)
        split_boundary = true
        smooth = true        
    elseif testCase == 2
        r = 2.0
        F,V = platonicsolid(4,r)     
        split_boundary = false
        smooth = true
    elseif testCase == 3 
        boxDim = [3.0, 3.0, 3.0] # Dimensions for the box in each direction
        pointSpacing = 1.0
        F,V,C = tribox(boxDim,pointSpacing)
        split_boundary = false
        smooth = true    
    elseif testCase == 4
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = [TriangleFace{Int}(f) for f in faces(M)]
        V = coordinates(M)
        F, V = mergevertices(F,V)
        split_boundary = false
        smooth = true
    end

    n = 1
    constrain_boundary = false
    F_tri, V_tri = subtri_dual(F, V, n; smooth=smooth, split_boundary=split_boundary, constrain_boundary=constrain_boundary)

    ## Visualization
    strokewidth1 = 0.5
    lineWidth = 4

    Fs, Vs = separate_vertices(F,V)

    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1], title = "Original")
    meshplot!(ax1, Fs, Vs, strokewidth=strokewidth1, strokecolor=:red, color=:white, transparency=false)

    # meshplot!(ax1, GeometryBasics.Mesh(V_dual,F_dual), strokewidth=strokewidth1, color=:red, transparency=true)

    ax2 = AxisGeom(fig[1, 2], title = "Refined, n=$n, smooth=$smooth, split_boundary=$split_boundary")
    hp = meshplot!(ax2, F_tri, V_tri, strokewidth=strokewidth1, strokecolor=:red, color=:white, transparency=false)

    # meshplot!(ax1, Fn1, Vn1, strokewidth=strokewidth1)
    # hp =scatter!(ax2, V_tri[end-length(V)-1:end], color=:red, markersize=15)

    n_range = 0:1:6
    hSlider2 = Slider(fig[2, :], range = n_range, startvalue = n, linewidth=30)
    on(hSlider2.value) do n    
        F_tri, V_tri = subtri_dual(F, V, n; smooth=smooth, split_boundary=split_boundary, constrain_boundary=constrain_boundary)
        hp[1] = GeometryBasics.Mesh(V_tri, F_tri)    
        ax2.title = "Refined, n=$n, smooth=$smooth, split_boundary=$split_boundary"
    end

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end