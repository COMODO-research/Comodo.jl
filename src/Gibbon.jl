module Gibbon

import GeometryBasics
import LinearAlgebra
import XML
import DataStructures
import FEBio
import Rotations

# import GLMakie

export gibbonDir, elements2indices, hexMeshBox, minDist, distND, lerp, gridPoints, interp_biharmonicSpline, interp_biharmonic_ND, gunique, ind2sub,sub2ind, meshEdges, unique_simplices, unique_dict, midEdgePoints, subtri, geoSphere, icosahedron, tetrahedron, cube, dodecahedron, octahedron, platonicsolid, meshnormal, toGeometryBasicsMesh, toGeometryBasicsSimplices, toGeometryBasicsPoints 

include("functions.jl")

end # module

