using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Rotations

GLMakie.closeall()

voxelSize = (2.0, 2.0, 2.0)

pointSpacing = 2.0 
rb = 60.0
nc = 150
V1 = [Point{3, Float64}(rb*cos(t), rb*sin(t), 0.0) for t in range(0.0, 1.3*pi, nc)]
R1 = collect(range(40.0, 20.0, nc))

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
nz = Point{3, Float64}(0.0, 0.0, 1.0)
P_cut_vec = normalizevector.([cross(V1[1], nz), cross(nz, V1[end])]) # Vectors pointing to distance origin
# P_cut_vec = [normalizevector(V1[1]-V1[2]), normalizevector(V1[end]-V1[end-1])] # Vectors pointing to distance origin
P_cut_origins = [V1[1], V1[end]] # Origins for vectors
D_cut_vec = [π.*R1[1], π.*R1[end]] # Distance from origin to consider for cut

FMc, VMc = cutends(FM, VM, P_cut_origins, P_cut_vec, D_cut_vec)

Eb = boundaryedges(FMc)
G = meshgroup(Eb)
numGroups = maximum(G)

## Visualization
fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Graph with radial data and derived surface")
hl = scatter!(ax1, V1, markersize = 10, color = R1, colorrange=(0.0, rMax), colormap=:viridis)
hm1 = meshplot!(ax1, FM, VM, color=(:white, 0.5), strokewidth=0.0, transparency= true)
Colorbar(fig[1, 2], hl)

scatter!(ax1, P_cut_origins, markersize = 25, color = :black)
arrows3d!(ax1, P_cut_origins, [R1[1], R1[end]].*P_cut_vec, color = :black)

ax2 = AxisGeom(fig[1, 3], title = "Cut surface")
hm1 = meshplot!(ax2, FMc, VMc, color=:white, strokewidth=0.1)

hm2 = edgeplot!(ax2, Eb, VMc, color=:blue, linewidth=3.0, depth_shift=-0.05f0)

screen = display(GLMakie.Screen(), fig)