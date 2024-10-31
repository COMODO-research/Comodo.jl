module Comodo

# Import dependancy libraries
import GeometryBasics
import LinearAlgebra
import DataStructures
import Rotations
import Statistics
import GLMakie
import Interpolations # E.g. for resampling curves
import BSplineKit
import QuadGK
import Distances
import DelaunayTriangulation # For regiontrimesh
import StaticArrays # For volumetric mesh definitions

include("functions.jl")

# Export imported packages
export GeometryBasics

# Export types
export Element, Tet4, Tet10, Hex8, Hex20, Penta6, Rhombicdodeca14, Truncatedocta24 # Volumetric elements (polyhedra)

# Export types/structs
export ConnectivitySet

# Export functions
export comododir, elements2indices, hexbox, mindist, dist, lerp
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
export circlepoints, loftlinear, trisurfslice
export edgeangles, count_edge_face, boundaryedges, boundaryfaces, boundaryfaceindices, edges2curve
export pointspacingmean, extrudecurve, meshgroup, ray_triangle_intersect
export distmarch, mesh_curvature_polynomial, curve_length, evenly_sample, evenly_space
export invert_faces, kabsch_rot, sweeploft, loftpoints2surf, revolvecurve
export batman, tridisc, quaddisc, regiontrimesh, scalesimplex, subcurve, dualclad
export tet2hex, element2faces, subhex, rhombicdodecahedron, tri2quad
export tetgenmesh, surfacevolume, tetvolume, extrudefaces, filletcurve
export squircle, circlerange, edgefaceangles, faceanglesegment, eulerchar
export rhombicdodecahedronfoam, kelvinfoam, truncatedoctahedron, ntrapezohedron, hexagonaltrapezohedron #, tetrakaidecahedron
export mag, indexmap!, indexmap, minp, maxp

# Export functions: Visualization related
export slidercontrol, slider2anim, dirplot, normalplot

end # module

