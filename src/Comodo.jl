module Comodo

# Import required functions and modules from dependency libraries
using Statistics: mean, Statistics
using Distances: euclidean, Distances
using QuadGK: quadgk, QuadGK
using StaticArrays: StaticVector, Size, StaticArrays
using Rotations: RotMatrix3, RotXYZ, rotation_between, AngleAxis, Rotations
using DataStructures: OrderedDict, DataStructures
import MarchingCubes # For isosurface creation
using TetGen: tetrahedralize, TetGen
using BSplineKit: BSplineOrder, BSplineKit
using DelaunayTriangulation: triangulate, each_solid_triangle, get_points, DelaunayTriangulation
using GLMakie: Slider, Axis3, Figure, LScene, Keyboard, events, record, set_close_to!, wireframe!, GLMakie
using LinearAlgebra: cross, norm, dot, eigen, svd, det, LinearAlgebra
using GeometryBasics: LineFace, Point, NgonFace, 
                      OffsetInteger, AbstractPoint, Vec, 
                      QuadFace, TriangleFace, faces, 
											coordinates, Vec3, GeometryBasics

include("functions.jl")

# Export imported modules for later possible use
export GeometryBasics
export Statistics
export Distances
export QuadGK
export StaticArrays
export Rotations
export MarchingCubes
export TetGen
export BSplineKit
export DelaunayTriangulation
export GLMakie
export LinearAlgebra
export GeometryBasics
export DataStructures

# Export (finite) element types
export AbstractElement, TetrahedronElement, PentahedronElement, HexahedronElement, TruncatedoctahedronElement, RhombicdodecahedronElement
export Tet4, Tet10, Tet15, Hex8, Hex20, Penta6, Penta15, Rhombicdodeca14, Truncatedocta24 

# Export types/structs
export ConnectivitySet

# Export functions
export comododir, slidercontrol, slider2anim, dirplot, normalplot,
elements2indices, hexbox, mindist, dist, lerp,
gridpoints, gridpoints_equilateral, 
interp_biharmonic_spline, interp_biharmonic, nbezier, 
ind2sub, sub2ind, meshedges, edgelengths, 
gunique, unique_simplices, unique_dict, occursonce, 
subtri, subquad, geosphere,hemisphere,quad2tri,  
icosahedron, tetrahedron, cube, dodecahedron, octahedron, platonicsolid, 
tofaces, topoints, togeometrybasics_mesh, 
triplate, quadplate, quadsphere, smoothmesh_laplacian, smoothmesh_hc, 
mergevertices, separate_vertices, remove_unused_vertices, 
edgecrossproduct, facearea, facenormal, vertexnormal, normalizevector, 
con_edge_face, con_edge_edge, con_face_edge, con_face_face, con_face_face_v, con_vertex_edge, con_vertex_vertex_f, 
con_vertex_face, con_vertex_vertex, con_vertex_simplex, meshconnectivity, 
simplexcenter, vertex2simplexdata, simplex2vertexdata, 
circlepoints, loftlinear, grid2surf, trisurfslice, 
edgeangles, count_edge_face, boundaryedges, boundaryfaces, boundaryfaceindices, edges2curve, 
pointspacingmean, extrudecurve, meshgroup, ray_triangle_intersect, 
distmarch, mesh_curvature_polynomial, curve_length, evenly_sample, evenly_space, 
invert_faces, invert_faces!, kabsch_rot, sweeploft, revolvecurve, 
batman, tridisc, quaddisc, regiontrimesh, scalesimplex, subcurve, dualclad, 
tet2hex, element2faces, subhex, rhombicdodecahedron, tri2quad, 
tetgenmesh, surfacevolume, tetvolume, extrudefaces, filletcurve, 
squircle, circlerange, edgefaceangles, faceanglesegment, eulerchar, 
rhombicdodecahedronfoam, kelvinfoam, truncatedoctahedron, ntrapezohedron, hexagonaltrapezohedron,  #, tetrakaidecahedron
mag, indexmap!, indexmap, minp, maxp, spacing2numvertices, 
joingeom, quadbox, tribox, tetbox, pad3, getisosurface, 
randangle, stepfunc, perlin_noise, removepoints, inpolygon, 
elementEdges, tet4_tet10, penta6_penta15, 
findindexin, hexagonline, hexagongrid, hexagonmesh, fromtomesh, fromtomesh!,
vectorpair_angle, triangulateboundary 

end # module