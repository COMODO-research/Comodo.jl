using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `bary2cart` to compute the Cartesian coordinates
`p` for the input barycentric point `λ` and the input triangle `f` and vertices `V`. 
=#

GLMakie.closeall()

function lambda2colorvec(λ)
    return convert(Vector{Float64}, λ)
end

for testCase = 1:2
     if testCase == 1
        # Create example triangle
        f, V = equilateraltriangle(2.0) # Single equilateral triangle domain

        # Bary centric coordinates
        λ = Point{3,Float64}(1/3,1/3,1/3)

        # Compute Cartesian coordinates
        p = bary2cart(f,V,λ)

        # Visualisation ----------------------------------------------------------------
        function updateBaryPlot(λ₁, λ₂, λ₃)
            λ = Point{3,Float64}(λ₁, λ₂, λ₃)
            hp2.color = lambda2colorvec(λ)
            ax.title = "λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3]) * ", Σ=" *string(sum(λ))
            hp4[1] = bary2cart(f, V, λ)
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

        stepRange = -1.0:0.001:1.0
        hSlider1 = Slider(fig[2, :], range = stepRange, startvalue = 1/3, linewidth=30)

        on(hSlider1.value) do λ₁ # Get 1 from slider  
            λ₂ = hSlider2.value[] 
            λ₃ = hSlider3.value[]    
            updateBaryPlot(λ₁, λ₂, λ₃)
        end

        hSlider2 = Slider(fig[3, :], range = stepRange, startvalue = 1/3, linewidth=30)
        on(hSlider2.value) do λ₂     
            λ₁ = hSlider1.value[] 
            λ₃ = hSlider3.value[]    
            updateBaryPlot(λ₁, λ₂, λ₃)
        end

        hSlider3 = Slider(fig[4, :], range = stepRange, startvalue = 1/3, linewidth=30)
        on(hSlider3.value) do λ₃   
            λ₁ = hSlider1.value[] 
            λ₂ = hSlider2.value[]
            updateBaryPlot(λ₁, λ₂, λ₃)
        end

        Colorbar(fig[:, 2], hp2)

        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")

    elseif testCase == 2
        e = Tet4{Int64}(1, 2, 3, 4)
          
        r = 1.0 

        # Create vertices    
        a = r*sqrt(2.0)/sqrt(3.0)
        b = -r*sqrt(2.0)/3.0
        c = -r/3.0       

        V=Vector{Point{3, Float64}}(undef,4)
        V[1 ] = Point{3, Float64}(   -a,      b,   c)
        V[2 ] = Point{3, Float64}(    a,      b,   c)    
        V[3 ] = Point{3, Float64}(  0.0,    0.0,   r)
        V[4 ] = Point{3, Float64}(  0.0, -2.0*b,   c)

        f = element2faces(e)

        # Bary centric coordinates
        λ = Point{4,Float64}(1/3, 1/3, 1/3, 1/3)

        # Compute Cartesian coordinates
        p = bary2cart(e, V, λ)

        # Visualisation ----------------------------------------------------------------
        function updateBaryPlot3(λ₁, λ₂, λ₃, λ₄)
            λ = Point{4,Float64}(λ₁, λ₂, λ₃, λ₄) 
            hp2.color = lambda2colorvec(λ)
            ax.title = "λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3]) * ", λ₄="*string(λ[4]) * ", Σ=" *string(sum(λ))
            hp4[1] = bary2cart(e, V, λ)
        end

        cmap = :bluesreds # Colormap to use for bary centric coordinates

        c = convert(Vector{Float64},λ)
        fig2 = Figure(size = (800,800))
        ax = AxisGeom(fig2[1, 1], title="λ₁="*string(λ[1])*", λ₂="*string(λ[2])*", λ₃="*string(λ[3])*", λ₄="*string(λ[4]))
        hp1 = meshplot!(ax, f, V; color = (:gray,0.1), strokewidth=3, transparency=true)
        hp2 = scatter!(ax, V, markersize=35, color=c, depth_shift = -0.01f0, colormap=cmap, colorrange=(-1.0,1.0))
        hp3 = text!(ax, V; text = ["P1","P2","P3", "P4"], font = :bold, fontsize=25, depth_shift=-0.1f0)
        hp4 = scatter!(ax, p, markersize=25, color=:black, depth_shift = -0.01f0)


        stepRange = -1.0:0.01:1.0

        tolLevel = 1e-2
        hSlider1 = Slider(fig2[2, :], range = stepRange, startvalue = 1/3, linewidth=30)
        
        on(hSlider1.value) do λ₁ # Get 1 from slider  
            λ₂ = hSlider2.value[] 
            λ₃ = hSlider3.value[]
            λ₄ = hSlider4.value[]
            updateBaryPlot3(λ₁, λ₂, λ₃, λ₄)
        end

        hSlider2 = Slider(fig2[3, :], range = stepRange, startvalue = 1/3, linewidth=30)
        on(hSlider2.value) do λ₂ # Get 1 from slider  
            λ₁ = hSlider1.value[] 
            λ₃ = hSlider3.value[]
            λ₄ = hSlider4.value[]
            updateBaryPlot3(λ₁, λ₂, λ₃, λ₄)
        end

        hSlider3 = Slider(fig2[4, :], range = stepRange, startvalue = 1/3, linewidth=30)
        on(hSlider3.value) do λ₃ # Get 1 from slider  
            λ₁ = hSlider1.value[]
            λ₂ = hSlider2.value[]             
            λ₄ = hSlider4.value[]
            updateBaryPlot3(λ₁, λ₂, λ₃, λ₄)
        end

        hSlider4 = Slider(fig2[5, :], range = stepRange, startvalue = 1/3, linewidth=30)
        on(hSlider4.value) do λ₄ # Get 1 from slider  
            λ₁ = hSlider1.value[]
            λ₂ = hSlider2.value[] 
            λ₃ = hSlider3.value[]
            updateBaryPlot3(λ₁, λ₂, λ₃, λ₄)
        end

        Colorbar(fig2[:, 2], hp2)

        screen = display(GLMakie.Screen(), fig2)
        GLMakie.set_title!(screen, "testCase = $testCase")
    end

end