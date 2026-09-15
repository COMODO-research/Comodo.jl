using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using FileIO
using Printf

GLMakie.closeall()

for testCase  in 1:6
    if testCase == 1
        r = 10.0
        F, V = geosphere(5,r) 
        z = [v[3] for v ∈ V]
        indStart = findmax(z)[2]   
        indEnd = findmin(z)[2] # indEnd = rand(1:lastindex(V),1)[1]  
        con_V2V = con_vertex_vertex_f(F, V) 
    elseif testCase == 2
        r = 1.0
        n = 3 # Number of refinement steps from cube
        F, V, C = subquadsphere(n, r)
        z = [v[3] for v ∈ V]
        indStart = findmax(z)[2]   
        indEnd = findmin(z)[2] # indEnd = rand(1:lastindex(V),1)[1]  
        con_V2V = con_vertex_vertex_f(F, V) 
    elseif testCase == 3
        plateDim = [10.0, 10.0]
        plateElem = [10, 10]
        orientation = :up
        F, V, Eb, Cb = quadplate(plateDim, plateElem; orientation=orientation)
        dn = [norm(v) for v ∈ V]
        indStart = findmax(dn)[2]   
        indEnd = findmin(dn)[2] # indEnd = rand(1:lastindex(V),1)[1]  
        E = meshedges(F)
        con_V2V = con_vertex_vertex(E, V) 
    elseif testCase == 4
        plateDim = [10.0, 10.0]
        plateElem = [10, 10]
        orientation = :up
        F, V, Eb, Cb = quadplate(plateDim, plateElem; orientation=orientation)
        dn = [norm(v) for v ∈ V]
        indStart = findmax(dn)[2]   
        indEnd = findmin(dn)[2] # indEnd = rand(1:lastindex(V),1)[1]          
        con_V2V = con_vertex_vertex_f(F, V) 
    elseif testCase == 5   
        w = 3.0     
        boxDim = [w, w, w] # Dimensions for the box in each direction
        boxEl = [5, 5, 5] # Number of elements to use in each direction 
        F, V, C = quadbox(boxDim,boxEl)
        _, indStart = findmin(norm.(V.-Point{3, Float64}(boxEl[1]/2.0, boxEl[2]/2.0, boxEl[3]/2.0)))
        _, indEnd = findmin(norm.(V.-Point{3, Float64}(-boxEl[1]/2.0, -boxEl[2]/2.0, -boxEl[3]/2.0)))
        con_V2V = con_vertex_vertex_f(F, V) 
    elseif testCase == 6
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F, V, _, _ = mergevertices(F, V)
        z = [v[3] for v ∈ V]
        indStart = findmax(z)[2]   
        d, dd, _ = distmarch(F, V, [indStart]) # Compute distances marched from start 
        _, indEnd = findmax(d)
        con_V2V = con_vertex_vertex_f(F, V) 
    end   

    pathVec, distVec = meshgeodesic(F, V, indStart, indEnd; con_V2V=con_V2V)

    ## Visualization
    cmap = Makie.Reverse(:Spectral)
    d, _, _ = distmarch(F, V, [indStart]; con_V2V=con_V2V) # Compute distances marched from start 

    function nanmaximum(d)
        return maximum(x->isnan(x) ? -Inf : x, d)
    end

    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Distances")
    hp1 = meshplot!(ax1, F, V, color=d, strokewidth=0.1, colormap=cmap, colorrange=(0.0, nanmaximum(d)), transparency=true)
    
    # scatter!(ax1, V, color=d, markersize=15, depth_shift=-0.01f0, colormap=cmap, colorrange=(0.0, nanmaximum(d)))

    scatter!(ax1,V[indStart], color=:black, markersize=25, depth_shift=-0.01f0)
    scatter!(ax1,V[indEnd], color=:black, markersize=25, depth_shift=-0.01f0)

    # text!(ax1, V; text=[@sprintf("%.3f",d) for d in d], depth_shift=-0.02f0)

    # for k in keys(dd)
    #     println(k)
    #     text!(ax1, (0.5*V[k[1]] + 0.5* V[k[2]]); text=@sprintf("%.3f", dd[k]) , depth_shift=-0.02f0, color=:green)
    # end

    hp2 = lines!(ax1, V[pathVec], color=distVec, linewidth=5.0, depth_shift=-0.1f0, colormap=cmap, colorrange=(0.0, nanmaximum(d)))
    Colorbar(fig[1, 2], hp1)

    
    screen = display(GLMakie.Screen(), fig)

end