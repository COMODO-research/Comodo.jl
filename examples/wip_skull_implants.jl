using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.GLMakie.Colors

using Geogram
using FileIO

#=
This demo is for the `kabsch_rot` function. The Kabsch method allows one to 
determine the rotation that occurred between two meshes (or point sets) with 
point-to-point correspondence. In this demo a triangulated surface is loaded 
from an STL file. Next the coordinates are rotated, and the Kabsch method is 
used to retrieve the rotation performed, allowing one to "unrotate" the data.  
=#

pointSpacing = 5.0

cBone = RGBf(0.89, 0.855, 0.788)
cBoneA = RGBA(0.89, 0.855, 0.788, 0.25)
cPurple = RGBf(1.0, 0.30196078431372547, 0.023529411764705882)

# Loading a mesh
pathName = "/home/kevin/DATA/Code/Julia/BodyParts3D/assets/BodyParts3D_data/stl/"
fileNames = ["FMA52734.stl", # Frontal bone
            "FMA52735.stl", # Occipital bone
            "FMA52736.stl", # Sphenoid bone
            "FMA52738.stl", # Right temporal bone
            "FMA52739.stl", # Left temporal bone
            "FMA52788.stl", # Right parietal bone
            "FMA52789.stl", # Left parietal bone
            "FMA52892.stl", # Right zygomatic bone
            "FMA52893.stl", # Left zygomatic bone
            "FMA53645.stl", # Right lacrimal bone
            "FMA53646.stl", # Left lacrimal bone
            "FMA53647.stl", # Right nasal bone
            "FMA53648.stl", # Left nasal bone
            "FMA53649.stl", # Right maxilla
            "FMA53650.stl", # Left maxilla
            "FMA53655.stl", # Right palatine bone
            "FMA53656.stl", # Left palatine bone
            "FMA54737.stl", # Right inferior nasal concha
            "FMA54738.stl", # Left inferior nasal concha
            "FMA52748.stl", # Mandible
            "FMA55680.stl", # Right upper lateral incisor tooth
            "FMA55681.stl", # Right upper central incisor tooth
            "FMA55682.stl", # Left upper central incisor tooth
            "FMA55683.stl", # Left upper lateral incisor tooth
            "FMA55686.stl", # Right lower canine tooth
            "FMA55687.stl", # Left lower canine tooth
            "FMA55688.stl", # Right second upper premolar tooth
            "FMA55689.stl", # Right first upper premolar tooth
            "FMA55690.stl", # Left first upper premolar tooth
            "FMA55691.stl", # Left second upper premolar tooth
            "FMA55692.stl", # Left second lower premolar tooth
            "FMA55693.stl", # Left first lower premolar tooth
            "FMA55694.stl", # Right first lower premolar tooth
            "FMA55695.stl", # Right second lower premolar tooth
            "FMA55697.stl", # Right second upper molar tooth
            "FMA55698.stl", # Right first upper molar tooth
            "FMA55699.stl", # Left first upper molar tooth
            "FMA55700.stl", # Left second upper molar tooth 
            "FMA55703.stl", # Left second lower molar tooth
            "FMA55704.stl", # Left first lower molar tooth
            "FMA55705.stl", # Right first lower molar tooth
            "FMA55706.stl", # Right second lower molar tooth
            "FMA55798.stl", # Right upper canine tooth
            "FMA55799.stl", # Left upper canine tooth
            "FMA57140.stl", # Right lower lateral incisor tooth
            "FMA57141.stl", # Left lower lateral incisor tooth
            "FMA57142.stl", # Right lower central incisor tooth
            "FMA57143.stl", # Left lower central incisor tooth
             ] #joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")

function getModels(pathName::String, fileNames::Vector{String})
    FM = Vector{Vector{TriangleFace{Int}}}(undef, length(fileNames))
    VM = Vector{Vector{Point{3, Float64}}}(undef, length(fileNames))
    for (i, fileName) in enumerate(fileNames)
        M = load(joinpath(pathName, fileName))
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F, V, _, _ = mergevertices(F, V)
        FM[i] = F
        VM[i] = V
    end
    F, V, C = joingeom(FM, VM)
    return F, V, C 
end

F, V, C = getModels(pathName, fileNames)

## Define bone to replace 

replaceLabel = 8

function getbone(F, V, C, replaceLabel)
    # Get bone to replace
    B_rep = C.==replaceLabel
    F_rep = F[B_rep]
    F_rep, V_rep = remove_unused_vertices(F_rep, V)
    return F_rep, V_rep, B_rep
end

F_rep, V_rep, B_rep = getbone(F, V, C, replaceLabel)

# Remesh using geomgram
numPoints = spacing2numvertices(F_rep, V_rep, pointSpacing)
F_repr, V_repr = ggremesh(F_rep, V_rep; nb_pts=numPoints)#, remesh_anisotropy=1.0, remesh_gradation = 1.0, pre_max_hole_area=100, pre_max_hole_edges=0, post_max_hole_area=100, post_max_hole_edges=0, quiet=0, suppress = true)


## Visualization
GLMakie.closeall()
cmap = Makie.Categorical(:Spectral)

Fs, Vs = separate_vertices(F, V)
Cs = simplex2vertexdata(Fs, C, Vs)

F_repr_s, V_repr_s = separate_vertices(F_repr, V_repr)

fig = Figure(size = (1400,1000))

ax1 = AxisGeom(fig[1, 1])
hp1 = meshplot!(ax1, Fs, Vs; color=Cs, strokewidth=0.1, colormap=cmap)
Colorbar(fig[1, 2], hp1)

ax2 = AxisGeom(fig[1, 3])
hp2 = meshplot!(ax2, Fs[.!B_rep], Vs; color=cBoneA, strokewidth=0.0, transparency=true)
hp3 = meshplot!(ax2, F_repr_s, V_repr_s; color=(:grey, 1.0), strokewidth=0.25, strokecolor=:gray)

screen = display(GLMakie.Screen(), fig)

## 

function dualclad_solid(F, V, s)
    Fs, Fq, Vs = dualclad(F,V,s; connectivity=:face)

    Ns = vertexnormal(Fs, Vs)

    t = pointSpacing/5.0
    E1, V1 = extrudefaces(Fs, Vs; extent=t, direction=:positive, num_steps=2, N=Ns);
    F_E1 = element2faces(E1)
    E2, V2 = extrudefaces(Fq, Vs; extent=t, direction=:positive, num_steps=2, N=Ns);
    E2, V2 = subhex(E2, V2, 1; direction=2)
    F_E2 = element2faces(E2)

    numVertices_V1 = length(V1)
    V_dual = deepcopy(V1)
    append!(V_dual, V2)
    Fq = [QuadFace{Int}(f.+numVertices_V1) for f in F_E2]
    F_dual_quad = F_E1[2]
    append!(F_dual_quad, Fq)
    F_dual_quad, V_dual = mergevertices(F_dual_quad, V_dual)
    indBoundaryFaces = boundaryfaceindices(F_dual_quad)
    F_dual_quad = F_dual_quad[indBoundaryFaces]
    F_dual_quad, V_dual = remove_unused_vertices(F_dual_quad, V_dual)

    F_dual_quad, V_dual = subquad(F_dual_quad, V_dual, 1; method=:linear)
    F1q, V1q = tri2quad(F_E1[1], V1)
    numVertices_V_dual = length(V_dual)
    append!(F_dual_quad, [QuadFace{Int}(f.+numVertices_V_dual) for f in F1q])
    append!(V_dual, V1q)
    F_dual_quad, V_dual = mergevertices(F_dual_quad, V_dual)
    F_dual_quad, V_dual = subquad(F_dual_quad, V_dual, 2; method=:Catmull_Clark)

    n=25
    λ=0.5
    V_dual = smoothmesh_laplacian(F_dual_quad, V_dual, n, λ)

    return F_dual_quad, V_dual
end

s = 0.3

F_dual_quad, V_dual = dualclad_solid(F_repr, V_repr, s)

##

fig2 = Figure(size = (1400,1000))

ax1 = AxisGeom(fig2[1, 1])
hp1 = meshplot!(ax1, Fs[.!B_rep], Vs; color=cBone, strokewidth=0.0)#, transparency=true)
hp2 = meshplot!(ax1, F_dual_quad, V_dual; color=cPurple, strokewidth=0.25, strokecolor=:gray)


stepRange1 = range(0.1, 0.9, 25)
hSlider1 = Slider(fig2[2, :], range = stepRange1, startvalue = s, linewidth=30)

on(hSlider1.value) do s
    replaceLabel = hSlider2.value[]
    F_rep, V_rep, B_rep = getbone(F, V, C, replaceLabel)

    # Remesh using geomgram
    numPoints = spacing2numvertices(F_rep, V_rep, pointSpacing)
    F_repr, V_repr = ggremesh(F_rep, V_rep; nb_pts=numPoints)#,

    F_dual_quad, V_dual = dualclad_solid(F_repr, V_repr, s)

    hp2[1] = GeometryBasics.Mesh(V_dual, F_dual_quad)
end

stepRange2 = 1:20#length(fileNames)
hSlider2 = Slider(fig2[3, :], range = stepRange2, startvalue = replaceLabel, linewidth=30)

on(hSlider2.value) do replaceLabel

    F_rep, V_rep, B_rep = getbone(F, V, C, replaceLabel)

    # Remesh using geomgram
    numPoints = spacing2numvertices(F_rep, V_rep, pointSpacing)
    F_repr, V_repr = ggremesh(F_rep, V_rep; nb_pts=numPoints)#, remesh_anisotropy=1.0, remesh_gradation = 1.0, pre_max_hole_area=100, pre_max_hole_edges=0, post_max_hole_area=100, post_max_hole_edges=0, quiet=0, suppress = true)


    s = hSlider1.value[]

    F_dual_quad, V_dual = dualclad_solid(F_repr, V_repr, s)

    hp1[1] = GeometryBasics.Mesh(Vs, Fs[.!B_rep])
    hp2[1] = GeometryBasics.Mesh(V_dual, F_dual_quad)
end


screen = display(GLMakie.Screen(), fig2)
