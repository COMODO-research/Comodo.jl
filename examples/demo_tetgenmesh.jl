using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using FileIO

testCase = 6

if testCase == 1
    r = 1.0
    Fb,Vb = geosphere(3,r)
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb)    
elseif testCase == 2
    r1 = 2.0
    r2 = r1/2
    F1,V1 = geosphere(3,r1)
    F2,V2 = geosphere(2,r2)
    F2 = [f.+length(V1) for f in invert_faces(F2)]    
    Fb = [F1;F2]
    Vb = [V1;V2]

    D = edgelengths(Fb,Vb) 
    vol1 = mean(D)^3 / (6.0*sqrt(2.0))

    Cb = ones(length(Fb))
    Cb[length(F1)+1:end] .+= 1

    v_hole = Point{3,Float64}(0.0, 0.0, 0.0)
    v_region = Point{3,Float64}(r2+((r1-r2)/2), 0.0, 0.0)

    stringOpt = "paAqY"
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=[v_region],region_vol=vol1,V_holes=[v_hole], stringOpt)
elseif testCase == 3
    r1 = 2.0
    r2 = r1/2
    F1,V1 = geosphere(3,r1)
    D1 = edgelengths(F1,V1) 

    F2,V2 = geosphere(3,r2)
    D2 = edgelengths(F2,V2) 

    F2 = [f.+length(V1) for f in invert_faces(F2)]    
    Fb = [F1;F2]
    Vb = [V1;V2]
    
    vol1 = mean(D1)^3 / (6.0*sqrt(2.0))
    vol2 = mean(D2)^3 / (6.0*sqrt(2.0))

    Cb = ones(length(Fb))
    Cb[length(F1)+1:end] .+= 1

    v_region1 = Point{3,Float64}(r2+((r1-r2)/2), 0.0, 0.0)
    v_region2 = Point{3,Float64}(0.0, 0.0, 0.0)    

    region_vol = [vol1,vol2]

    stringOpt = "paAqY"
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=[v_region1,v_region2],region_vol=region_vol, stringOpt)
elseif testCase == 4
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)
    F1 = tofaces(faces(M))
    V1 = [Point{3,Float64}(v) for v in topoints(coordinates(M))]            
    F1,V1,ind1,ind2 = mergevertices(F1,V1; roundVertices=true)

    v_hole = Point{3,Float64}(10.0,-15.0,-25.0)

    F2,V2 = geosphere(3,30)    
    V2 .*= Point{3,Float64}(1.5,1.0,1.0)
    V2 .+= v_hole
    F2 = [f.+length(V1) for f in invert_faces(F2)]

    Fb = [F1;F2]
    Vb = [V1;V2]
    
    D = edgelengths(Fb,Vb) 
    vol1 = mean(D)^3 / (6.0*sqrt(2.0))

    Cb = ones(length(Fb))
    Cb[length(F1)+1:end] .+= 1

    v_region = Point{3,Float64}(-50.0, 0.0, 0.0)

    stringOpt = "paAqYQ"
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=[v_region],V_holes=[v_hole], region_vol=vol1, stringOpt)
elseif testCase == 5
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)
    F1 = tofaces(faces(M))
    V1 = [Point{3,Float64}(v) for v in topoints(coordinates(M))]            
    F1,V1,ind1,ind2 = mergevertices(F1,V1; roundVertices=true)    
    vol1 = mean(edgelengths(F1,V1))^3 / (6.0*sqrt(2.0))

    v_cent = Point{3,Float64}(10.0,-15.0,-25.0)

    F2,V2 = geosphere(4,30)    
    V2 .*= Point{3,Float64}(1.5,1.0,1.0)
    V2 .+= v_cent
    vol2 = mean(edgelengths(F2,V2))^3 / (6.0*sqrt(2.0))
    
    F2 = [f.+length(V1) for f in F2]
    
    F3,V3 = geosphere(3,20)        
    V3 .+= v_cent
    F3 = [f.+length(V1).+length(V2) for f in invert_faces(F3)]
    

    Fb = [F1;F2;F3]
    Vb = [V1;V2;V3]

    Cb = ones(length(Fb))
    Cb[length(F1)+1:end] .+= 1
    Cb[length(F1)+length(F2)+1:end] .+= 1

    v_region1 = Point{3,Float64}(-50.0, 0.0, 0.0)
    v_region2 = Point{3,Float64}(10.0,25-15.0,-25.0)    
    V_regions = [v_region1, v_region2]

    V_holes = [v_cent]

    region_vol = [vol1,vol2]

    stringOpt = "paAqYQ"
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=V_regions,V_holes=V_holes, region_vol=region_vol, stringOpt)
elseif testCase == 6
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)
    F1 = tofaces(faces(M))
    V1 = [Point{3,Float64}(v) for v in topoints(coordinates(M))]            
    F1,V1,ind1,ind2 = mergevertices(F1,V1; roundVertices=true)    
    vol1 = mean(edgelengths(F1,V1))^3 / (6.0*sqrt(2.0))

    v_cent = Point{3,Float64}(10.0,-15.0,-25.0)

    F2,V2 = geosphere(4,30)    
    V2 .*= Point{3,Float64}(1.5,1.0,1.0)
    V2 .+= v_cent
    vol2 = mean(edgelengths(F2,V2))^3 / (6.0*sqrt(2.0))
    
    F2 = [f.+length(V1) for f in F2]
    
    F3,V3 = geosphere(3,20)        
    V3 .+= v_cent
    F3 = [f.+length(V1).+length(V2) for f in invert_faces(F3)]
    

    Fb = [F1;F2;F3]
    Vb = [V1;V2;V3]

    Cb = ones(length(Fb))
    Cb[length(F1)+1:end] .+= 1
    Cb[length(F1)+length(F2)+1:end] .+= 1

    v_region1 = Point{3,Float64}(-50.0, 0.0, 0.0)
    v_region2 = Point{3,Float64}(10.0,25-15.0,-25.0)    
    V_regions = [v_region1, v_region2]

    V_holes = [v_cent]

    region_vol = [vol1,vol2]

    stringOpt = "paAqYQ"
    element_type = Tet10{Int}
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb; facetmarkerlist=Cb, V_regions=V_regions,V_holes=V_holes, region_vol=region_vol, stringOpt, element_type=element_type)
end

## Visualization
cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # Triangular faces
CE_F = repeat(CE,inner=4)

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)

Fs,Vs = separate_vertices(F,V)
CE_Vs = simplex2vertexdata(Fs,CE_F)
M = GeometryBasics.Mesh(Vs,Fs)

strokewidth = 1 

fig = Figure(size=(1500,1500))

ax1 = Axis3(fig[1, 1][1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary surfaces")
hp1 = mesh!(ax1,GeometryBasics.Mesh(Vbs,Fbs), color=Cbs_V, shading = FastShading, transparency=true, overdraw=false,colorrange = (1,3),colormap=cmap)

scatter!(ax1,V_regions,color=:black,markersize=25)
scatter!(ax1,V_holes,color=:black,markersize=25)

ax2 = Axis3(fig[1, 1][1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Cut mesh")
hp2 = poly!(ax2,M, color=CE_Vs, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=strokewidth, overdraw=false,colorrange = (1,2),colormap=cmap)
hp3 = scatter!(ax2,V,color=:black,markersize=10)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F,V)))

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = mean(stepRange),linewidth=30)

on(hSlider.value) do z 
    B = ZE .<= z
    indShow = findall(B)
    if isempty(indShow)
        hp2.visible=false        
    else        
        hp2.visible=true
        Fs = element2faces(E[indShow])
        Cs = repeat(CE[indShow],inner=4)
        
        indB = boundaryfaceindices(Fs)        
        Fs = Fs[indB]
        Cs = Cs[indB]
        Fs,Vs = separate_vertices(Fs,V)
        CE_Vs = simplex2vertexdata(Fs,Cs)
        Ms = GeometryBasics.Mesh(Vs,Fs)
        hp2[1] = Ms
        hp2.color = CE_Vs
        hp3[1] = V[elements2indices(E[indShow])]
    end
end
# hSlider.selected_index[]+=1
slidercontrol(hSlider,ax2)

fig

# fileName = comododir()*"/assets/img/TetGen_example_01.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=6)

# fig




