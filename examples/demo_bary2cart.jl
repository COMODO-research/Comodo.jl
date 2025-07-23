using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `bary2cart` to compute the Cartesian coordinates
`p` for the input barycentric point `λ` and the input triangle `f` and vertices `V`. 
=#

GLMakie.closeall()

testCase = 1 

# Create example triangle
f, V = equilateraltriangle(2.0) # Single equilateral triangle domain

# Bary centric coordinates
λ = Point{3,Float64}(1/3,1/3,1/3)

# Compute Cartesian coordinates
p = bary2cart(f,V,λ)

# Visualisation ----------------------------------------------------------------
function lambda2colorvec(λ)
    return convert(Vector{Float64},λ)
end

function updateBaryPlot(λ₁,λ₂,λ₃)
    λ = Point{3,Float64}(λ₁,λ₂,λ₃)
    hp2.color = lambda2colorvec(λ)
    ax.title = "λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3]) * ", Σ=" *string(sum(λ))
    hp4[1] = bary2cart(f,V,λ)
end

cmap = :bluesreds # Colormap to use for bary centric coordinates

fig = Figure(size = (800,800))
ax = AxisGeom(fig[1, 1], azimuth=-pi/2, elevation=pi/2, 
title="λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3]), 
limits=(-2.0,2.0,-2.0,2.0,-2.0,2.0))

hp1 = meshplot!(ax, f, V)
hp2 = scatter!(ax, V, markersize=35, color=lambda2colorvec(λ), depth_shift = -0.01f0, colormap=cmap, colorrange=(-1.0,1.0))
hp3 = text!(ax, V; text = ["P1","P2","P3"], font = :bold, fontsize=25, depth_shift=-0.02f0)
hp4 = scatter!(ax, p, markersize=25, color=:red, depth_shift = -0.02f0)

stepRange = -1/3:1/12:1.0
hSlider1 = Slider(fig[2, :], range = stepRange, startvalue = 1/3, linewidth=30)

on(hSlider1.value) do λ₁    
    λ₂ = hSlider2.value[]
    λ₃ = 1.0 - λ₁ - λ₂ 
    if λ₃<minimum(stepRange)
        λ₃ = minimum(stepRange)
        λ₂ = 1.0 - λ₁ - λ₃
        set_close_to!(hSlider2, λ₂)
    end
    set_close_to!(hSlider3, λ₃)         
    updateBaryPlot(λ₁,λ₂,λ₃)
end

hSlider2 = Slider(fig[3, :], range = stepRange, startvalue = 1/3, linewidth=30)
on(hSlider2.value) do λ₂     
    λ₁ = hSlider1.value[]
    λ₃ = 1.0 - λ₁ - λ₂ 
    if λ₃<minimum(stepRange)
        λ₃ = minimum(stepRange)
        λ₁ = 1.0 - λ₂ - λ₃
        set_close_to!(hSlider1, λ₁)
    end
    set_close_to!(hSlider3, λ₃) 
    updateBaryPlot(λ₁,λ₂,λ₃)λ₂
end

hSlider3 = Slider(fig[4, :], range = stepRange, startvalue = 1/3, linewidth=30)
on(hSlider3.value) do λ₃   
    λ₂ = hSlider2.value[]
    λ₁ = 1.0 - λ₂ - λ₃
    if λ₁<minimum(stepRange)
        λ₁ = minimum(stepRange)
        λ₂ = 1.0 - λ₁ - λ₃
        set_close_to!(hSlider2, λ₂)
    end
    set_close_to!(hSlider1, λ₁) 
    updateBaryPlot(λ₁,λ₂,λ₃)
end

Colorbar(fig[:, 2], hp2)

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "testCase = $testCase")