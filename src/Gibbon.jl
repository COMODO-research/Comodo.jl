module Gibbon

# Import dependancy libraries
import GeometryBasics
import LinearAlgebra
import XML
import DataStructures
import FEBio
import Rotations
import Statistics
import SparseArrays
# import GLMakie

export indIn, gibbonDir, sliderControl, elements2indices, hexMeshBox, minDist, distND, lerp, gridPoints, interp_biharmonicSpline, interp_biharmonic_ND, gunique, ind2sub,sub2ind, meshEdges, unique_simplices, unique_dict, midPoints, subTri, subQuad, geoSphere, icosahedron, tetrahedron, cube, dodecahedron, octahedron, platonicsolid, meshnormal, toGeometryBasicsMesh, toGeometryBasicsSimplices, toGeometryBasicsPoints 
export mergeVertices, con_edge_face, con_edge_edge, con_face_edge, con_face_face, con_vertex_edge, con_vertex_face, con_vertex_vertex, meshCon, connectivitySet

include("functions.jl")

end # module

