module Comodo

# Import required functions and modules from dependency libraries
using Statistics: mean, Statistics
using Distances: euclidean, Distances, Euclidean, pairwise
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
export comododir, slidercontrol, slider2anim, dirplot, normalplot
export elements2indices, hexbox, mindist, dist, lerp
export gridpoints, gridpoints_equilateral
export interp_biharmonic_spline, interp_biharmonic, nbezier
export ind2sub, sub2ind, meshedges, edgelengths
export gunique, unique_simplices, unique_dict, occursonce
export subtri, subquad, geosphere,hemisphere,quad2tri
export icosahedron, tetrahedron, cube, dodecahedron, octahedron, platonicsolid
export tofaces, topoints, togeometrybasics_mesh
export triplate, quadplate, quadsphere, smoothmesh_laplacian, smoothmesh_hc
export mergevertices, separate_vertices, remove_unused_vertices
export edgecrossproduct, facearea, facenormal, vertexnormal, normalizevector
export con_edge_face, con_edge_edge, con_face_edge, con_face_face, con_face_face_v, con_vertex_edge, con_vertex_vertex_f
export con_vertex_face, con_vertex_vertex, con_vertex_simplex, meshconnectivity
export simplexcenter, vertex2simplexdata, simplex2vertexdata
export circlepoints, loftlinear, grid2surf, trisurfslice
export edgeangles, count_edge_face, boundaryedges, boundaryfaces, boundaryfaceindices, edges2curve
export pointspacingmean, extrudecurve, meshgroup, ray_triangle_intersect
export distmarch, mesh_curvature_polynomial, curve_length, evenly_sample, evenly_space
export invert_faces, invert_faces!, kabsch_rot, sweeploft, revolvecurve
export batman, tridisc, quaddisc, regiontrimesh, scalesimplex, subcurve, dualclad
export tet2hex, element2faces, subhex, rhombicdodecahedron, tri2quad
export tetgenmesh, surfacevolume, tetvolume, extrudefaces, filletcurve
export squircle, circlerange, edgefaceangles, faceanglesegment, eulerchar
export rhombicdodecahedronfoam, kelvinfoam, truncatedoctahedron, ntrapezohedron, hexagonaltrapezohedron  #, tetrakaidecahedron
export mag, indexmap!, indexmap, minp, maxp, spacing2numvertices
export joingeom, quadbox, tribox, tetbox, pad3, getisosurface
export randangle, stepfunc, perlin_noise, removepoints, inpolygon
export elementEdges, tet4_tet10, penta6_penta15
export findindexin, hexagonline, hexagongrid, hexagonmesh, fromtomesh, fromtomesh!
export vectorpair_angle, triangulateboundary, faceinteriorpoint

end # module