using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.Rotations
using Comodo.LinearAlgebra

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

###########################################################

pointSpacing = 2.0
r = 16.0
h = 70.0 # Height = extrusion distance (extent)

rBranch = r/2
hBranch = 50.0
filletDistance = rBranch/4

# Create base branch
nc = ceil(Int,(2*pi*r)/pointSpacing) # Compute number of circumferential points from point spacing
Vc = circlepoints(r,nc;dir=:acw) # Create points on circle

n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0)) # Extrusion direction
direction = :positive
F1,V1 = extrudecurve(Vc; extent=h, direction=:positive, n=n, close_loop=true,face_type=:quad)
E1b  = boundaryedges(F1)

F1t = quad2tri(F1,V1; convert_method = :angle);
N1_V = vertexnormal(F1,V1)
E1 = boundaryedges(F1)

# Define side branch vector and origin
tBranch = 0.25*π
Q = RotXYZ(tBranch,0.0,0.0)
nz = Vec{3,Float64}(0.0,0.0,1.0)
m = Q*nz

m = Vec{3, Float64}(0.0,sqrt(2)/2,sqrt(2)/2)
pb = Point{3,Float64}(0.0,0.0,h/4)

# Define side branch base contour 
nc_branch = ceil(Int,(2*pi*rBranch)/(pointSpacing)) # Compute number of circumferential points from point spacing
Vc_branch_c = circlepoints(rBranch,nc_branch;dir=:acw) # Create points on circle
# Q = rotation_between(m,nz)

b = 1.0 / cos(pi-tBranch);
Qb = RotXYZ(0.5*π,0.0,0.0)


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

β = [acosd(dot(m,n))/2.0 for n in N1_V[indNearest]] 
dCut = filletDistance ./ tand.(β)

Vc_branch_ellipse=[pb+(Qb*Point{3,Float64}(v[1], b*v[2], v[3]))+m*(r*(1/cos(tBranch))+maximum(dCut)) for v in Vc_branch_c];

# March distances from intersection point  
V1[indNearest] = P_intersect # Overwrite nearest with intersection points so distances are more accurate 
d,dd,l = distmarch(F1,V1,indNearest)

# Cut hole in surface 
boolKeep = d .> maximum(dCut)

F1_cut,F1_cutout = removefaces(F1,boolKeep)

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
Vc_bez = [pb+(Qb*Point{3,Float64}(v[1], b*v[2], v[3]))+m*(r*(1/cos(tBranch))+maximum(dCut)) for v in Vc_bez];

Vc_bez2 = circlepoints(rBranch,nc;dir=:cw)
Vc_bez2 = [pb+(Q'*v)+m*hBranch for v in Vc_bez2];

# Vc_bez2 = circlepoints(rBranch,nc;dir=:cw)
# Vc_bez2 = [pb+(Qb*Point{3,Float64}(v[1], b*v[2], v[3]))+m*hBranch for v in Vc_bez2];


_,iClosest = findmin(norm.(Vc_bez .- V1[indCurve[1]]))
if iClosest > 1 
    circshift!(Vc_bez,nc-iClosest+1)
end

numPointBezier = 8
v1 = filletDistance*4 # Departure velocity of Hermite equivalent spline
v2 = filletDistance*4 # Arrival velocity of Hermite equivalent spline
VV_B = []
Vb = Vector{Point{3,Float64}}()
nEnd = -m
for i in eachindex(indCurve)
    nBase = KV_cutout[i]    
    pb = [V1[indCurve[i]], V1[indCurve[i]] + v1/3*nBase, Vc_bez[i] + v2/3*nEnd,Vc_bez[i]]
    vb = nbezier(pb,numPointBezier)
    push!(VV_B,vb)
    append!(Vb,vb)
end

Fb = grid2surf(Vb,length(indCurve);periodicity=(false,true),face_type=:quad)
Eb = boundaryedges(Fb)
indB = reduce(vcat,Eb)
nSmooth = 25
α=0.1
β=0.5
Vb = smoothmesh_hc(Fb,Vb, nSmooth, α, β; constrained_points=indB)

Fb2,Vb2 = loftlinear(Vc_bez2,Vc_bez;num_steps=25,close_loop=true,face_type=:quad)

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

E,VE = extrudefaces(F,V; extent=2.0, direction=:negative, num_steps=3) 

FE = element2faces(E)

## Visualization
linewidth = 6
strokewidth = 0.5
markersize = 20
cmap = cgrad(:Spectral, 3, categorical = true)

fig = Figure(size=(1900,1200))
# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
# # ax1 = LScene(fig[1,1])

# hp1 = lines!(ax1,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
# hp2 = poly!(ax1,GeometryBasics.Mesh(V1,F1_cut), strokewidth=strokewidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap = :Spectral)
# hp3 = dirplot(ax1,pb,m,scaleval=hBranch,color=:black,linewidth=8)

# hp5 = lines!(ax1,Vc_bez,color=:red,linewidth=linewidth)
# hp5 = lines!(ax1,Vc_bez2,color=:blue,linewidth=linewidth)

# # hp5 = scatter!(ax1,P_intersect,color=:red,markersize=markersize)
# # hp6 = scatter!(ax1,V1[indNearest],color=:cyan,markersize=markersize)

# # dirplot(ax1,V1[indCurve],MV_cutout[indCurve],scaleval=2,color=:blue,linewidth=2)
# # dirplot(ax1,V1[indCurve],KV_cutout,scaleval=2,color=:red,linewidth=2)

# # hp6 = scatter!(ax1,V1[indNear],color=:blue,markersize=markersize)

# # hp6 = scatter!(ax1,V1[indMiddle],color=:cyan,markersize=markersize)


# # scatter!(ax1,Vc_branch,color=:red,markersize=markersize)

# wireframe!(ax1,GeometryBasics.Mesh(V1,E1_cutout),linewidth=linewidth, transparency=false, color=:yellow)

# hp2 = poly!(ax1,GeometryBasics.Mesh(Vb,Fb), strokewidth=strokewidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap = :Spectral)
# hp2 = poly!(ax1,GeometryBasics.Mesh(Vb2,Fb2), strokewidth=strokewidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false,colormap = :Spectral)

# # Colorbar(fig[1, 2], hp2)

# ax2 = LScene(fig[1,1])

Fs,Vs = separate_vertices(F,V)
Cs = simplex2vertexdata(Fs,C,Vs)

ax1 = AxisGeom(fig[1, 1])
hp11 = meshplot!(ax1,GeometryBasics.Mesh(Vs,Fs), strokewidth=strokewidth,color=Cs, strokecolor=:black, colormap = cmap,colorrange=(0.5,3.5), stroke_depth_shift=-0.001f0)
Colorbar(fig[1, 2], hp11)

FEs,VEs = separate_vertices(FE,VE)

ax2 = AxisGeom(fig[1, 3])
hp21 = meshplot!(ax2,GeometryBasics.Mesh(VEs,FEs), strokewidth=strokewidth,color=:white, strokecolor=:black, colormap = cmap,colorrange=(0.5,3.5), stroke_depth_shift=-0.001f0)

fig


# save("/home/kevin/DATA/Julia/Comodo.jl/assets/temp/branches.png", fig, px_per_unit=5)