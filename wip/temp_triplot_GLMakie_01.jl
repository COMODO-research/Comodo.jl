using GLMakie
using DelaunayTriangulation
using Random


# Random.seed!(1234)
# points = randn(Point2f, 50)
# tri = triangulate(points)


outer = [
    (0.0,0.0),(2.0,1.0),(4.0,0.0),
    (6.0,2.0),(2.0,3.0),(3.0,4.0),
    (6.0,6.0),(0.0,6.0),(0.0,0.0)
]
inner = [
    (1.0,5.0),(2.0,4.0),(1.01,1.01),
    (1.0,1.0),(0.99,1.01),(1.0,5.0)
]
boundary_points = [[outer], [inner]]
boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
tri = triangulate(points; boundary_nodes = boundary_nodes)
refine!(tri; max_area=1e-2*get_total_area(tri))


############# 
GLMakie.activate!(inline=false)
fig = Figure(fontsize=24)

ax = Axis(fig[1, 1], title = "(a): Unconstrained", titlealign = :left, width=400,height=400)
triplot!(ax, tri) 
fig