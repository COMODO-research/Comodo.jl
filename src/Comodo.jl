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

include("functions.jl")

# Export functions
export comododir, elements2indices, hexbox, mindist, dist, lerp
export gridpoints, interp_biharmonic_spline, interp_biharmonic, nbezier, gunique
export ind2sub, sub2ind, meshedges, edgelengths, unique_simplices, unique_dict
export subtri, subquad, geosphere, quad2tri 
export icosahedron, tetrahedron, cube, dodecahedron, octahedron, platonicsolid
export tofaces, topoints, togeometrybasics_mesh
export quadplate, quadsphere, smoothmesh_laplacian, smoothmesh_hc
export mergevertices, separate_vertices, remove_unused_vertices
export edgecrossproduct, facearea, facenormal, vertexnormal, normalizevector, dirplot, normalplot
export con_edge_face, con_edge_edge, con_face_edge, con_face_face, con_face_face_v, con_vertex_edge, con_vertex_vertex_f
export con_vertex_face, con_vertex_vertex, con_vertex_simplex, meshconnectivity
export simplexcenter, vertex2simplexdata, simplex2vertexdata
export circlepoints, loftlinear, trisurfslice
export wrapindex, edgeangles, count_edge_face, boundaryedges, edges2curve
export pointspacingmean, extrudecurve, meshgroup, ray_triangle_intersect
export distmarch, mesh_curvature_polynomial, curve_length, evenly_sample
export invert_faces, kabsch_rot, sweeploft, loftpoints2surf, revolvecurve
export batman

# Export functions: Visualization related
export slidercontrol

# Export types/structs
export ConnectivitySet

# Export imported packages
export GeometryBasics

end # module

