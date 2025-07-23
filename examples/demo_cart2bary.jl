using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `cart2bary` to compute the barycentric coordinates
`λ` for the input point `p` and the input triangle `f` and vertices `V`. 
=#

GLMakie.closeall()

testCase = 1
# Create example triangle

# Here a refined set of triangles is created and the example triangle is 
# located in the middle. This way we can use the other triangles to create a set 
# of query points for the visualisation below.  
f, V = equilateraltriangle(1.0) # Single equilateral triangle domain
Fn,Vn = subtri([f],V,1) # Split once into 4 triangles 
f = convert(TriangleFace{Int},Fn[1]) # Keep first (middle) triangle
f,V,_ = remove_unused_vertices(f,Vn) # Create clean point set of just the triangle

# Define point to compute barycentric coordinates for
p = Point{3,Float64}(0.0, 0.0, 0.0)

# Compute bary centric coordinates
λ = cart2bary(f,V,p)

# Visualisation ----------------------------------------------------------------
# This visualisation creates a set of points in and around the input triangle 
# so the user can inspect the barycentric coordinates for these by manipulating 
# slider 

# Define visualisation helper function that will create the triangles to 
# visualise the barycentric coordinate magnitudes onto. 
function createTriangleSet(f,V,p,λ)
     Vp = deepcopy(V)
     push!(Vp,p)
     n = length(Vp)
     Fp = [TriangleFace{Int}(f[2], f[3], n), 
           TriangleFace{Int}(f[3], f[1], n), 
           TriangleFace{Int}(f[1], f[2], n)]
     Fps,Vps = separate_vertices(Fp, Vp)
     Cps = simplex2vertexdata(Fps, convert(Vector{Float64},λ))     
     return Fps, Vps, Cps
end

Fps, Vps, Cps = createTriangleSet(f,V,p,λ)
_,PG = subtri(Fn,Vn,2) # Define point set as thrice refined 
pushfirst!(PG,Point{3,Float64}(0.0, 0.0, 0.0)) # Also add origin
c = convert(Vector{Float64},λ) # Color to use for visualization of scatter and triangles

cmap = :bluesreds # Colormap to use for bary centric coordinates

# Start visualisation
fig = Figure(size = (800,800))
ax = AxisGeom(fig[1, 1], azimuth=-pi/2, elevation=pi/2, title="λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3]))
hp1 = edgeplot!(ax, f, V; color = :black, linewidth=3, depth_shift = -0.02f0)
hp2 = meshplot!(ax, Fps, Vps; color = Cps, depth_shift = -0.01f0, stroke_depth_shift = -0.01f0, colormap=cmap, colorrange=(-1.0,1.0))
hp3 = scatter!(ax, p, markersize=25, color=:black, depth_shift = -0.01f0)
scatter!(ax, PG, markersize=10, color=:black, depth_shift = -0.02f0)
hp4 = scatter!(ax, V, markersize=35, color=c, depth_shift = -0.01f0, colormap=cmap, colorrange=(-1.0,1.0))

text!(ax, V; text = ["P1","P2","P3"], font = :bold, fontsize=25, depth_shift=-0.1f0)

Colorbar(fig[1, 2], hp2)

stepRange = 1:length(PG)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 1, linewidth=30)

on(hSlider.value) do i   
     p = PG[i]
     λ = cart2bary(f,V,p)
     Fps, Vps, Cps = createTriangleSet(f,V,p,λ)
     c = convert(Vector{Float64},λ)
     hp2[1] = GeometryBasics.Mesh(Vps,Fps)
     hp2.color = Cps
     hp3[1] = p 
     hp4.color = c
     ax.title = "λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3])
end
slidercontrol(hSlider,ax)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "testCase = $testCase")