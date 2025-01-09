using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.BSplineKit
using Comodo.QuadGK: quadgk # For numerical integration

testCase = 2

# Define input curve
if testCase == 1
    V = Point{3,Float64}[[0,0,0],[1.0,0,0],[1,1,0],[0,1,0]]
    close_loop = false
    c = 3
elseif testCase == 2
    n = 250
    r = 2*pi
    V = [Point{3,Float64}(r*cos(t),r*sin(t),0.0) for t in range(0,2*pi,n)]
    close_loop = true
    c = 2*pi*r
end

L = curve_length(V; close_loop = close_loop)
d = last(L)

println("Max. length: " * string(d) * ", theoretical/true: " * string(c))
