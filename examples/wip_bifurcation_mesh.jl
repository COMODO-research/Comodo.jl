using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.Rotations
using Comodo.LinearAlgebra

GLMakie.closeall()

function removefaces(F1,boolKeep)
    F1_cut = Vector{eltype(F1)}()
    F1_cutout = Vector{eltype(F1)}()
    for f in F1
        if all(boolKeep[f])
            push!(F1_cut,f)
        else
            push!(F1_cutout,f)
        end
    end
    return F1_cut,F1_cutout
end

function loopsmooth(V; λ=0.5)    
    nv = length(V)
    wSelf = (1.0-λ)
    wOther = λ/2.0
    for i in eachindex(V)
        V[i] = wOther*V[mod1(i-1,nv)] + wSelf*V[i] + wOther*V[mod1(i+1,nv)] 
    end
    return V
end

tBranch = 90.0
wallThickness = 8.0
r = 125.0
h = 500.0 # Height = extrusion distance (extent)
rBranch = 57.0
hBranch = 400.0
filletDistance = 10#rBranch/4
v1 = filletDistance*3 # Departure velocity of Hermite equivalent spline
v2 = v1 # Arrival velocity of Hermite equivalent spline

pointSpacing = 10.0 

# Create base branch quad mesh
nc = ceil(Int,(2*pi*r)/pointSpacing) # Compute number of circumferential points from point spacing
Vc = circlepoints(r,nc;dir=:acw) # Create points on circle
n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0)) # Extrusion direction
F1,V1 = extrudecurve(Vc; extent=h, direction=:positive, n=n, close_loop=true,face_type=:quad)

function addBranch(F,V, rBranch, hBranch, tBranch, filletDistance, v1, v2; numPointsBezierInitial = 100)

    F1 = deepcopy(F)
    V1 = deepcopy(V)

    pointSpacing = pointspacingmean(F1,V1)

    E1b  = boundaryedges(F1)

    # Convert quad mesh to triangulation for ray tracing
    F1t = quad2tri(F1,V1; convert_method = :angle);
    N1_V = vertexnormal(F1,V1)
    E1 = boundaryedges(F1)

    # Define side branch vector and origin
    Q = RotXYZ(tBranch*(π/180.0),0.0,0.0)
    nz = Vec{3,Float64}(0.0,0.0,1.0)
    m = Q'*nz

    # Define branch base point on main centre line
    pb = Point{3,Float64}(0.0,0.0, h/2.0-r/tand(tBranch))

    # Define side branch base contour 
    nc_branch = ceil(Int,(2*pi*rBranch)/(pointSpacing)) # Compute number of circumferential points from point spacing
    Vc_branch_c = circlepoints(rBranch,nc_branch;dir=:acw) # Create points on circle
    # Q = rotation_between(m,nz)

    b = 1.0 / sind(180.0-tBranch) # Ellipse scaling factor
    Qb = RotXYZ(-0.5*π,0.0,0.0)

    Vc_branch = [pb+(Q'*v)+m*hBranch for v in Vc_branch_c]

    # Project branch contour onto main
    P_intersect = Vector{Point{3,Float64}}(undef,nc_branch)
    indIntersect = Vector{Int}(undef,nc_branch)
    boolShot = fill(false,length(F1t))
    for i in 1:nc_branch
        p,j = ray_triangle_intersect(F1t,V1,Vc_branch[i],-m; rayType = :ray, triSide = 1)
        P_intersect[i] = p[1]
        indIntersect[i] = j[1]
        boolShot[j[1]] = true
    end
    indIntersect_vertices = reduce(vcat,F1t[indIntersect])
    pc,jc = ray_triangle_intersect(F1t,V1,pb,m; rayType = :ray, triSide = -1)

    _,indMin = mindist(P_intersect,V1[indIntersect_vertices]; getIndex=Val(true))
    indNearest = indIntersect_vertices[indMin]

    β = [90.0-acosd(dot(m,n)) for n in N1_V[indNearest]] # Angle between branch direction vector and local normal
    dCut = filletDistance ./ sind.(β/2.0)
    dCutMax = maximum(dCut)
    Vc_branch_ellipse=[pb+(Qb*Point{3,Float64}(v[1], b*v[2], v[3]))+m*(r/sind(tBranch)+dCutMax) for v in Vc_branch_c];

    # March distances from intersection point  
    V1[indNearest] = P_intersect # Overwrite nearest with intersection points so distances are more accurate 
    d,dd,l = distmarch(F1,V1,indNearest)

    # Cut hole in surface 
    boolKeep = d .> dCutMax

    F1_cut, F1_cutout = removefaces(F1,boolKeep)

    C = meshgroup(F1_cut; con_type = :v, indStart=E1b[1][1], stop_at = 1)
    F1_cut = F1_cut[C.==1]

    E1b_cut = boundaryedges(F1_cut)
    bool_E1b_cutout = [!in(e,E1b) for e in E1b_cut]

    E1_cutout = E1b_cut[bool_E1b_cutout]
    indCurve = edges2curve(E1_cutout;remove_last = true)

    ME1_cutout = [normalizevector(V1[e[2]]-V1[e[1]]) for e in E1_cutout]
    MV_cutout  = simplex2vertexdata(E1_cutout,ME1_cutout,V1)

    VE1_cutout = vertex2simplexdata(E1_cutout,V1)
    NE1_cutout = vertex2simplexdata(E1_cutout,N1_V)

    KV_cutout = cross.(MV_cutout[indCurve],N1_V[indCurve])


    V1[indCurve] .= loopsmooth(V1[indCurve]; λ=1/3)    

    KV_cutout .= normalizevector(loopsmooth(KV_cutout; λ=0.5))
    KV_cutout .= normalizevector(loopsmooth(KV_cutout; λ=0.5))

    nc = length(indCurve)
    Vc_bez = circlepoints(rBranch,nc;dir=:cw)
    Vc_bez = [pb+(Qb*Point{3,Float64}(v[1], b*v[2], v[3]))+m*(r/sind(tBranch)+dCutMax) for v in Vc_bez];

    Vc_bez2 = circlepoints(rBranch,nc;dir=:cw)
    Vc_bez2 = [pb+(Q'*v)+m*hBranch for v in Vc_bez2];

    _,iClosest = findmin(norm.(Vc_bez .- V1[indCurve[1]]))
    if iClosest > 1 
        circshift!(Vc_bez,nc-iClosest+1)
    end

    # Create initial Bezier splines (not sample with correct spacing yet)
    VV_B = Vector{Vector{Point{3,Float64}}}()
    bezierLengths = Vector{Float64}(undef,length(indCurve))
    for i in eachindex(indCurve)
        nBase = KV_cutout[i]    
        pb = [V1[indCurve[i]], V1[indCurve[i]] + v1/3*nBase, Vc_bez[i] - v2/3*m,Vc_bez[i]]    
        vb = nbezier(pb, numPointsBezierInitial)    
        bezierLengths[i] = curve_length(vb; close_loop = false, total=true)
        push!(VV_B,vb)    
    end

    # Resample curves evenly based on base curve point spacing
    pointSpacingStart = pointspacingmean(V1[indCurve]; close_loop=true)
    pointSpacingEnd = pointspacingmean(Vc_bez; close_loop=true)
    numPointsBezier = ceil(Int,maximum(bezierLengths)/pointSpacingStart)
    Vb = Vector{Point{3,Float64}}()
    for i in eachindex(indCurve)
        vb = evenly_sample(VV_B[i], numPointsBezier; spline_order=4, close_loop = false, niter = 10)
        append!(Vb,vb)
    end

    Fb = grid2surf(Vb,length(indCurve); periodicity=(false,true),face_type=:quad)
    Eb = boundaryedges(Fb)
    indB = reduce(vcat,Eb)
    nSmooth = 25
    α=0.1
    β=0.5
    Vb = smoothmesh_hc(Fb,Vb, nSmooth, α, β; constrained_points=indB)

    num_steps_branch_straight = ceil(Int,minimum([norm(Vc_bez2[i]-Vc_bez[i]) for i in eachindex(Vc_bez)])/pointSpacingEnd)
    Fb2,Vb2 = loftlinear(Vc_bez2,Vc_bez; num_steps=num_steps_branch_straight,close_loop=true,face_type=:quad)

    F,V,C = joingeom(F1_cut,V1,Fb,Vb,Fb2,Vb2)
    F,V = remove_unused_vertices(F,V)
    F,V = mergevertices(F,V)

    ind1 = reduce(vcat,F[C.==1])
    ind2 = reduce(vcat,F[C.==2])
    ind3 = reduce(vcat,F[C.==3])

    con_F2F = con_face_face_v(F)

    boolSmooth = C.==2 #fill(true,length(F))
    for q=1:2
        for i in findall(boolSmooth)    
            boolSmooth[con_F2F[i]] .= true
        end
    end

    indConstrain = reduce(vcat,F[.!boolSmooth])

    nSmooth = 25
    α=0.1
    β=0.5
    V = smoothmesh_hc(F,V, nSmooth, α, β; constrained_points=indConstrain)

    return F, V, C
end

function addBranchHex(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2, wallThickness; numPointsBezierInitial = 100)

    F, V, C = addBranch(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2; numPointsBezierInitial = numPointsBezierInitial)

    # Extrude faces to hex8 elements 
    num_steps = 1+ceil(Int, wallThickness/pointSpacing)
    E, Ve = extrudefaces(F,V; extent=wallThickness, direction=:negative, num_steps=num_steps)
    elementLabels = repeat(C, outer=num_steps-1)
    elementLabels_F = repeat(elementLabels,inner=6) 
    FE = element2faces(E)
    indBoundaryFaces = boundaryfaceindices(FE)
    Fb = FE[indBoundaryFaces]
    Cb = elementLabels_F[indBoundaryFaces]
    return E, Ve, FE, Fb, Cb, F, V, C
end

E, Ve, FE, Fb, Cb, F, V, C = addBranchHex(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2, wallThickness; numPointsBezierInitial = 100)

Fbs,Vbs = separate_vertices(Fb,Ve)
Cbs = simplex2vertexdata(Fbs,Cb,Vbs)

## Visualization
linewidth = 6
strokewidth = 0.5
markersize = 20
cmap = cgrad(:Spectral, 3, categorical = true)

fig = Figure(size=(1900,1200))

Fs,Vs = separate_vertices(F,V)
Cs = simplex2vertexdata(Fs,C,Vs)

ax1 = AxisGeom(fig[1, 1])
hp1 = meshplot!(ax1, Fbs, Vbs, strokewidth=strokewidth, color=Cbs, strokecolor=:black, colormap = cmap, stroke_depth_shift=-0.001f0)
Colorbar(fig[1, 2], hp11)

stepRange1 = range(30.0,150.0, 100)
hSlider1 = Slider(fig[2, :], range = stepRange1, startvalue = tBranch,linewidth=30)
on(hSlider1.value) do tBranch 
    v1 = hSlider2.value[]
    filletDistance = hSlider3.value[]
    v2 = v1
    wallThickness = hSlider4.value[]
    E, Ve, FE, Fb, Cb, F, V, C = addBranchHex(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2, wallThickness; numPointsBezierInitial = 100)
    Fbs,Vbs = separate_vertices(Fb,Ve)
    Cbs = simplex2vertexdata(Fbs,Cb,Vbs)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)
    hp1.color = Cbs
end

stepRange2 = range(filletDistance/10,filletDistance*10, 500)
hSlider2 = Slider(fig[3, :], range = stepRange2, startvalue = v1,linewidth=30)
on(hSlider2.value) do v1
    tBranch = hSlider1.value[]
    filletDistance = hSlider3.value[]
    v2 = v1
    wallThickness = hSlider4.value[]
    E, Ve, FE, Fb, Cb, F, V, C = addBranchHex(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2, wallThickness; numPointsBezierInitial = 100)
    Fbs,Vbs = separate_vertices(Fb,Ve)
    Cbs = simplex2vertexdata(Fbs,Cb,Vbs)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)
    hp1.color = Cbs
end

stepRange3 = range(filletDistance/10,filletDistance*3, 500)
hSlider3 = Slider(fig[4, :], range = stepRange2, startvalue = filletDistance,linewidth=30)
on(hSlider3.value) do filletDistance
    tBranch = hSlider1.value[]
    v1 = hSlider2.value[]
    v2 = v1
    wallThickness = hSlider4.value[]
    E, Ve, FE, Fb, Cb, F, V, C = addBranchHex(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2, wallThickness; numPointsBezierInitial = 100)
    Fbs,Vbs = separate_vertices(Fb,Ve)
    Cbs = simplex2vertexdata(Fbs,Cb,Vbs)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)
    hp1.color = Cbs
end

stepRange4 = range(wallThickness/10,wallThickness*3, 100)
hSlider4 = Slider(fig[5, :], range = stepRange2, startvalue = wallThickness, linewidth=30)
on(hSlider4.value) do wallThickness
    tBranch = hSlider1.value[]
    v1 = hSlider2.value[]
    v2 = v1
    filletDistance = hSlider3.value[]    
    E, Ve, FE, Fb, Cb, F, V, C = addBranchHex(F1,V1, rBranch, hBranch, tBranch, filletDistance, v1, v2, wallThickness; numPointsBezierInitial = 100)
    Fbs,Vbs = separate_vertices(Fb,Ve)
    Cbs = simplex2vertexdata(Fbs,Cb,Vbs)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)
    hp1.color = Cbs
end


screen = display(GLMakie.Screen(), fig)


# ## Visualization
# linewidth = 6
# strokewidth = 0.5
# markersize = 20
# cmap = cgrad(:Spectral, 3, categorical = true)

# fig = Figure(size=(1900,1200))

# Fbs,Vbs = separate_vertices(Fb,Ve)
# Cbs = simplex2vertexdata(Fbs,Cb,Vbs)

# ax1 = AxisGeom(fig[1, 1])
# hp11 = meshplot!(ax1, Fbs, Vbs, strokewidth=strokewidth, color=Cbs, strokecolor=:black, colormap = cmap, stroke_depth_shift=-0.001f0)
# Colorbar(fig[1, 2], hp11)

# screen = display(GLMakie.Screen(), fig)

