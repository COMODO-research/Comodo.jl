"""
This demo shows the use of interp_biharmonic_spline for curve interpolation. 
"""

using Comodo
using GeometryBasics
using GLMakie
using LinearAlgebra
using BSplineKit
using QuadGK

# Define raw curve
testCase = 3
if testCase == 1
    np = 50
    t = range(0.0,2.0*π,np)
    V = [GeometryBasics.Point{3, Float64}(tt,3*sin(tt),0.0) for tt ∈ t]
elseif testCase == 2
    using Random
    Random.seed!(16)
    np = 50
    t = sort(unique(rand(np)*2.0*π)) 
    V = [GeometryBasics.Point{3, Float64}(tt,3*sin(tt),0.0) for tt ∈ t]
elseif testCase == 3
    np = 50
    t = range(0.0,2*π-(2.0*π/np),np)
    r = 2.0
    V = [GeometryBasics.Point{3, Float64}(r*cos(tt),sin(tt),sin(3.0*tt)) for tt ∈ t]
end

function evenly_sample(V,n)    
    m = length(V)
    T = curve_length(V) # Initialise as along curve (multi-linear) distance
    T ./= last(T) # Normalise
    S = interpolate(T, V, BSplineOrder(4),Natural()) # Create interpolator
    dS = Derivative() * S  # spline derivative
    L = zeros(Float64,m) # Initialise spline length vector
    for i ∈ 2:m
        # Compute length of segment [i-1, i]
        segment_length, _ = quadgk(T[i-1], T[i]; rtol = 1e-8) do t
            norm(dS(t))  # integrate |S'(t)| in segment [i, i + 1]
        end        
        L[i] = L[i - 1] + segment_length
    end    
    L ./= last(L) # Normalise to 0-1 range
    l = range(0.0,1.0,n) #Even allong curve distance 
    S = interpolate(L, V, BSplineOrder(4),Natural()) # Create interpolator
    return S.(l),S # Evaluate interpolator at even distance increments
end

# Evenly sample curve
n = 100; # Number of points for spline 
Vi,S = evenly_sample(V,n)

# Visualization
fig = Figure(size = (1200,800))

ax1 = Axis3(fig[1, 1],aspect = :data)
hp1 = scatter!(ax1, V,markersize=15,color=:red)
hp2 = lines!(ax1, V,linewidth=3,color=:red)

ax2 = Axis3(fig[1, 2],aspect = :data)
hp1 = scatter!(ax2, Vi,markersize=15,color=:black)
hp2 = lines!(ax2, Vi,linewidth=3,color=:black)
fig
