using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Rotations

GLMakie.closeall()

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

# Now computed radius weighted signed distance field
M = wsdf(xr, yr, zr, V1, R1; closest_type=:weighted)
siz = size(M)

# Compute level set surface for visualisation
FM, VM = getisosurface(M; x=xr, y=yr, z=zr, level=0.0, cap=false, padValue=1e8)      
# VM = smoothmesh_hc(FM, VM, 25)

## Cut branches

# Define cutting vectors and origins e.g. graph end points and end directions
n11 = length(V11)
n22 = length(V22)
n33 = length(V33)
n44 = length(V44)
indCut = [1, n11, n11+n22, n11+n22+n33, n11+n22+n33+n44]
v1 = normalizevector(V1[indCut[1]] -  V1[indCut[1]+1])
v2 = normalizevector(V1[indCut[2]] -  V1[indCut[2]-1])
v3 = normalizevector(V1[indCut[3]] -  V1[indCut[3]-1])
v4 = normalizevector(V1[indCut[4]] -  V1[indCut[4]-1])
v5 = normalizevector(V1[indCut[5]] -  V1[indCut[5]-1])

P_cut_vec = [v1, v2, v3, v4, v5] # Vectors pointing to distance origin
P_cut_origins = V1[indCut] # Origins for vectors
D_cut_vec = π.*R1[indCut] # Distance from origin to consider for cut

FMc, VMc = cutfrom(FM, VM, P_cut_origins, P_cut_vec, D_cut_vec)

## Visualization
fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Graph with radial data and derived surface")
hl = scatter!(ax1, V1, markersize = 10, color = R1, colorrange=(0.0, rMax), colormap=:viridis)
hm1 = meshplot!(ax1, FM, VM, color=(:white, 0.5), strokewidth=0.0, transparency= true)
Colorbar(fig[1, 2], hl)

scatter!(ax1, P_cut_origins, markersize = 25, color = :black)
arrows3d!(ax1, P_cut_origins, 30.0 .*P_cut_vec, color = :black)

ax2 = AxisGeom(fig[1, 3], title = "Cut surface")
hm1 = meshplot!(ax2, FMc, VMc, color=:white, strokewidth=0.1)

screen = display(GLMakie.Screen(), fig)