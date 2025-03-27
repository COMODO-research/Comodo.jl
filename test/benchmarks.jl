using Chairmarks
using Comodo
using Comodo.DelaunayTriangulation
Chairmarks.DEFAULTS.seconds = 5

function _elements2indices(F)
    return unique(reduce(vcat, F))
end

elements2indices_b = @b each_triangle(triangulate(rand(2, 10_000))) elements2indices,_elements2indices
# 1.328 ms (25 allocs: 192.625 KiB)
# 978.017 ms (59954 allocs: 4.474 GiB, 22.12% gc time)