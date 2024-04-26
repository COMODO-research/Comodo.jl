# Call required packages

using GeometryBasics # For point and mesh format
using LinearAlgebra # For things like dot and cross products
using DataStructures # For unique_dict
using Statistics # For: mean etc.
using GLMakie # For slidercontrol
using Rotations 
using Interpolations # E.g. for resampling curves
using BSplineKit # E.g. for resampling curves
using QuadGK: quadgk # For numerical integration
using Distances
using DelaunayTriangulation # For triangular meshing

"""
    ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V, con_V2V_f, con_F2F_v)

# Description 

A struct featuring the connectivity data for a mesh.   
"""
struct ConnectivitySet
    edge_vertex::Vector{LineFace{Int64}}
    edge_face::Vector{Vector{Int64}}
    edge_edge::Vector{Vector{Int64}}
    face_vertex#::Vector{Vector{Int64}} # Could be triangle/quad etc
    face_edge::Vector{Vector{Int64}}        
    face_face::Vector{Vector{Int64}}
    vertex_edge::Vector{Vector{Int64}}
    vertex_face::Vector{Vector{Int64}}        
    vertex_vertex::Vector{Vector{Int64}}
    vertex_vertex_f::Vector{Vector{Int64}}
    face_face_v::Vector{Vector{Int64}}
    ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V, con_V2V_f, con_F2F_v) = new(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V, con_V2V_f, con_F2F_v) 
end

"""
    comododir()

# Description 

This function simply returns the string for the Comodo path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function comododir()
    joinpath(@__DIR__, "..")
end


"""
    slidercontrol(hSlider,ax)    

Adds arrow key control to sliders

# Description 

This function adds arrow key control to GLMakie sliders. The inputs are the 
slider handle `hSlider` as well as the axis `ax`. If this function is called the
slider can be advanced a step by pressing the right arrow, and returned one step 
by pressing the left arrow. When one presses and holds the right or left arrow 
key, the slider will continue to move (as fast as graphics updating is possible
on your system) up to the end or start slider position respectively. Users may 
also use the up or down arrow keys. These function the same as the right and 
left arrow keys, however, rather than stopping at the slider extrema, the 
sliders position will "wrap" back to the start when advancing beyond the end 
position, and vice versa. 
"""
function slidercontrol(hSlider::Slider,ax::Union{Axis3, Figure})    
    sliderRange = hSlider.range[] # Get slider range
    rangeLength = length(sliderRange) # Number of possible steps 
    sliderIndex = hSlider.selected_index[] # Current slider index
    on(events(ax).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat # Pressed or held for instance            
            if event.key == Keyboard.up                                                  
                if sliderIndex == rangeLength
                    sliderIndex = 1                    
                else               
                    sliderIndex += 1                     
                end            
            elseif event.key == Keyboard.down                
                if sliderIndex == 1
                    sliderIndex = rangeLength
                else                
                    sliderIndex -= 1                    
                end                        
            end            
            if event.key == Keyboard.right                        
                if sliderIndex < rangeLength                       
                    sliderIndex += 1                     
                end            
            elseif event.key == Keyboard.left            
                if sliderIndex > 1                       
                    sliderIndex -= 1                     
                end            
            end
            hSlider.selected_index = sliderIndex            
        end
    end    
end

"""
    slider2anim(fig::Figure,hSlider::Slider,fileName::String; backforth=true, duration=2)

Exports movies from slider based visualisations.

# Description
Converts the effect of the slider defined by the slider handle `hSlider` for the 
figure `fig` to an animation/movie file
"""
function slider2anim(fig::Figure,hSlider::Slider,fileName::String; backforth=true, duration=2)
    stepRange = collect(hSlider.range[])
    if backforth
        append!(stepRange,reverse(stepRange[2:end-1]))
    end
    framerate = ceil(Int64,length(stepRange)/duration)
    record(fig, fileName, stepRange; framerate = framerate) do stepIndex
        set_close_to!(hSlider, stepIndex)
    end
end

"""
    elements2indices(F)

Returns the indices contained in `F`

# Description
    
This function obtains the unique set of indices for the vertices (nodes) 
used by the the simplices defined by `F`. The vector `F` may contain any type of 
simplices. For instance the elements in `F` may be of the type 
`GeometryBasics.TriangleFace` or `GeometryBasics.QuadFace` (or any other) for 
surface mesh data. However, volumetric elements of any type are permitted. In
essence this function simply returns `unique(reduce(vcat,F))`.
Hence any suitable vector containing vectors of numbers permitted by 
`reduce(vcat,F)` is supported. 
"""
function elements2indices(F)
    return unique(reduce(vcat,F))
end

"""
    gridpoints(x::Vector{T}, y=x, z=x) where T<:Real

Returns 3D grids of points

# Description

The `gridpoints` function returns a vector of 3D points which span a grid in 3D 
space. Points are defined as per the input ranges or range vectors. The output 
point vector contains elements of the type `Point`. 
"""
function gridpoints(x::Union{Vector{T}, AbstractRange{T}}, y=x, z=x) where T<:Real
    # reshape([Point{3, Float64}(x, y, z) for x in x, y in y, z in z],length(x)*length(y)*length(z)) # Features more allocations    
    V = Vector{Point{3,Float64}}(undef,length(x)*length(y)*length(z)) # Allocate point vector 
    c = 1 # Initiate linear index into V
    # Create grid with point order x->y->z (important for related meshing functions which assume this order)
    @inbounds for z in z  
        @inbounds for y in y
            @inbounds for x in x                           
                V[c] = Point{3,Float64}(x,y,z)
                c += 1 
            end
        end
    end
    return V
end  

"""
    gridpoints_equilateral(xSpan,ySpan,pointSpacing::T; return_faces = false, rectangular=false) where T <: Real

Returns a "grid" of 3D points that are located on the corners of an equilateral triangle tesselation.

# Description

This function returns 3D point data in the form of a `Vector{Point{3,Float64}}`. 
The point distribution is for an equilateral triangle tesselation. The input 
consists of the span in the x-, and y-direction, i.e. `xSpan` and `ySpan` 
respectively, as well as the desired `pointSpacing`. The "spans" should be 
vectors or tuples defining the minimum and maximum coordinates for the grid. The
true point spacing in the x-direction is computed such that a nearest whole 
number of steps can cover the required distance. Next this spacing is used to 
create the equilateral triangle point grid. Although the `xSpan` is closely 
adhered to through this method, the `ySpan` is not fully covered. In the 
y-direction the grid does start at the minimum level, but may stop short of 
reaching the maximum y as it may not be reachable in a whole number of steps 
from the minimum. Optional arguments include `return_faces` (default is 
`false`), which will cause the function to return triangular faces `F` as well 
as the vertices `V`. Secondly the option `rectangular` will force the grid to 
conform to a rectangular domain. This means the "jagged" sides are forced to be 
flat such that all x-coordinates on the left are at the minimum in `xSpan` and
all on the right are at the maximum in `xSpan`, however, this does result in a 
non-uniform spacing at these edges.  
"""
function gridpoints_equilateral(xSpan::Union{Vector{TT},Tuple{TT,TT}},ySpan::Union{Vector{TT},Tuple{TT,TT}},pointSpacing::T; return_faces = false, rectangular=false) where T <: Real where TT <: Real
    minX = minimum(xSpan)
    maxX = maximum(xSpan)
    minY = minimum(ySpan)
    maxY = maximum(ySpan)

    wx = maxX - minX
    nx = ceil(Int64,wx./pointSpacing)+1 # Number of points in the x-direction
    xRange = range(minX,maxX,nx) # The xRange
    pointSpacingReal_X=xRange[2]-xRange[1]

    pointSpacingReal_Y=pointSpacingReal_X.*0.5*sqrt(3)
    yRange = minY:pointSpacingReal_Y:maxY

    # Create the point grid
    V = Vector{Point{3,Float64}}(undef,length(xRange)*length(yRange))    
    c = 1
    sx = pointSpacingReal_X/4
    for j in eachindex(yRange)
        for i in eachindex(xRange)                 
            if iseven(j) # Shift over every second row of points
                x = xRange[i]+sx
            else
                x = xRange[i]-sx
            end                    
            if rectangular
                if isone(i)
                    x = minX
                elseif i == nx
                    x = maxX
                end
            end

            V[c] = Point{3}{Float64}(x,yRange[j],0.0)
            c += 1
        end
    end

    # Creat output, including faces if requested
    if return_faces == true
        plateElem=[length(xRange)-1,length(yRange)-1]
        F = Vector{TriangleFace{Int64}}(undef,prod(plateElem)*2)
        num_x = length(xRange)
        ij2ind(i,j) = i + ((j-1)*num_x) # function to convert subscript to linear indices    
        c = 1
        for i = 1:plateElem[1]
            for j = 1:plateElem[2]      
                if iseven(j)  
                    F[c] = TriangleFace{Int64}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1)])
                    c += 1
                    F[c] = TriangleFace{Int64}([ij2ind(i+1,j+1),ij2ind(i,j+1),ij2ind(i,j)])
                    c += 1
                else
                    F[c] = TriangleFace{Int64}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i,j+1)])
                    c += 1
                    F[c] = TriangleFace{Int64}([ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)])
                    c += 1
                end
            end
        end        
        return F,V
    else
        return V
    end
end

"""
    interp_biharmonic_spline(x::Union{Vector{T}, AbstractRange{T}},y::Union{Vector{T}, AbstractRange{T}},xi::Union{Vector{T}, AbstractRange{T}}; extrapolate_method=:linear,pad_data=:linear) where T<:Real

Interpolates 1D (curve) data using biharmonic spline interpolation

# Description

This function uses biharmonic spline interpolation [1], which features radial basis 
functions. The input is assumed to represent ordered data, i.e. consequtive 
unique points on a curve. The curve x-, and y-coordinates are provided through 
the input parameters `x` and `y` respectively. The third input `xi` defines the 
sites at which to interpolate. Each of in the input parameters can be either a 
vector or a range. 

# References
1. [David T. Sandwell, _Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data_, Geophysical Research Letters, 2, 139-142, 1987. doi: 10.1029/GL014i002p00139](https://doi.org/10.1029/GL014i002p00139)
"""
function interp_biharmonic_spline(x::Union{Vector{T}, AbstractRange{T}},y::Union{Vector{T}, AbstractRange{T}},xi::Union{Vector{T}, AbstractRange{T}}; extrapolate_method=:linear,pad_data=:linear) where T<:Real

    # Pad data if needed
    if isa(x,AbstractRange{T})
        xx = collect(x)
    else
        xx = deepcopy(x)
    end

    if isa(y,AbstractRange{T})
        yy = collect(y)
    else
        yy = deepcopy(y)
    end

    if pad_data==:linear
        # Linearly extended ends are added 
        dx1 = x[1]-x[2]
        dy1 = y[1]-y[2]    
        dx2 = x[end]-x[end-1]
        dy2 = y[end]-y[end-1]
           
        pushfirst!(xx,x[1]+dx1)
        pushfirst!(xx,x[1]+2.0*dx1) 
        push!(xx,x[end]+dx2)
        push!(xx,x[end]+2.0*dx2)
    
        pushfirst!(yy,y[1]+dy1)
        pushfirst!(yy,y[1]+2.0*dy1) 
        push!(yy,y[end]+dy2)
        push!(yy,y[end]+2.0*dy2)
    elseif pad_data==:constant
        # The start and end are copied twice 
        dx1 = x[1]-x[2]       
        dx2 = x[end]-x[end-1]
            
        pushfirst!(xx,x[1]+dx1)
        pushfirst!(xx,x[1]+2.0*dx1) 
        push!(xx,x[end]+dx2)
        push!(xx,x[end]+2.0*dx2)
    
        pushfirst!(yy,y[1])
        pushfirst!(yy,y[1]) 
        push!(yy,y[end])
        push!(yy,y[end])
    elseif  pad_data==:none
        # No padding
    else
        error("InvalidParameter: Invalid pad_data method provided, valued options are :linear, :constant, and :none")
    end

    # Change behaviour depending on extrapolation method
    if extrapolate_method==:linear
        # Simple data based linear extrapolation 
        L = [xii<x[1] || xii>x[end] for xii in xi] # Boolean for points to extrapolate for
        if any(L) # If any points outside of the range were encountered
            yi = Vector{Float64}(undef,length(xi)) # Initialise yi
            yi[L] = lerp(xx,yy,xi[L]) # Linearly extrapolate outside of the input range
            yi[.!L] = interp_biharmonic(xx,yy,xi[.!L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic(xx,yy,xi)
        end
    elseif extrapolate_method==:constant
        # Simple constant extrapolation (last/first value is repeated indefinately)
        Ls = xi.<x[1] # Boolean for points preceeding the start
        Ll = xi.>x[end] # Boolean for points after the end
        if any(Ls .|| Ll) # If any points outside of the range were encountered
            yi = Vector{Float64}(undef,length(xi)) # Initialise yi
            if any(Ls) # If any before start were found 
                yi[Ls] .= yy[1] # Just copy the start for these
            end
            if any(Ll) # If any after the end were found
                yi[Ll] .= yy[end] # Just copy the end for these
            end
            L = .!Ls .&& .!Ll # Boolean for points to interpolate using biharmonic interpolation
            yi[L] = interp_biharmonic(xx,yy,xi[L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic(xx,yy,xi)
        end
    elseif extrapolate_method==:biharmonic
        # Allow extrapolation as per the biharmonic function
        yi = interp_biharmonic(xx,yy,xi) 
    else
        error("InvalidParameter: Invalid extrapolate_method method provided, valued options are :linear, :constant, and :biharmonic")
    end

    return yi
end


"""
    interp_biharmonic(x,y,xi)

Interpolates n-dimensional data using biharmonic spline interpolation

# Description 

This function uses biharmonic interpolation [1]. The input `x` should define a 
vector consisting of m points which are n-dimensional, and the input `y` should
be a vector consisting of m scalar data values. 

# References 
1. [David T. Sandwell, _Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data_, Geophysical Research Letters, 2, 139-142, 1987. doi: 10.1029/GL014i002p00139](https://doi.org/10.1029/GL014i002p00139)
"""
function interp_biharmonic(x,y,xi)
    # Distances from all points in X to all points in X
    Dxx = dist(x,x)

    # Determine weights for interpolation
    g = (Dxx.^2) .* (log.(Dxx).-1.0)   # Green's function.
    g[1:size(g,1)+1:length(g)] .= 0.0 # Fix values along diagonal
    g[isnan.(g)] .= 0.0 # Replace NaN entries by zeros 
    W = g \ y # Weights  

    D = dist(xi,x) # Distance between points in X and XI

    G = (D.^2).*(log.(D).-1.0) # Green's function.
    G[isnan.(G)] .= 0.0 # Replace NaN entries by zeros
    return G * W
end

"""
    nbezier(P,n)

Returns a Bezier spline for the control points P whose order matches the numbe 
of control points provided. 

# Description 

This function returns `n` points for an m-th order Bézier spline, based on the 
m control points contained in the input vector `P`. This function supports point
vectors with elements of the type `AbstractPoint{3}` (e.g.
`Point{3, Float64}`) or `Vector{Float64}`.
"""
function nbezier(P::Vector{Point{ND,TV}},n::Integer) where ND where TV<:Real
    if n<2
        error("The vale of n is too low. Request at least two data points")
    end
    t = range(0,1,n) # t range
    N = length(P) # Number of control points 
    nn = 0:N-1
    nnr = N-1:-1:0
    f = factorial.(nn) 
    s = factorial(N-1)./( f.*reverse(f) ); # Sigma

    V =  zeros(Point{ND,TV},n)    
    @inbounds for i in 1:n
        b = s.* ((1.0-t[i]).^nnr) .* (t[i].^nn)
        @inbounds for j = 1:N
            V[i] += P[j].*b[j]
        end
    end 
    return V
end


"""
    lerp(x::Union{T,Vector{T}, AbstractRange{T}},y,xi::Union{T,Vector{T}, AbstractRange{T}}) where T <: Real

Linear interpolation

# Description 

This linearly interpolates (lerps) the input data specified by the sites `x` and 
data `y` at the specified sites `xi`. 
"""
function lerp(x::Union{T,Vector{T}, AbstractRange{T}},y,xi::Union{T,Vector{T}, AbstractRange{T}}) where T <: Real
    # Check if the lengths of x and y are the same
    if length(x)!=length(y)
        throw(DimensionMismatch("The number of sites in x should match the number of data entries in y"))
    end

    if length(xi)==1 # Single site provided
        yi = lerp_(x,y,xi)
    else # Loop over all data sites
        yi = Vector{eltype(y)}(undef,length(xi))
        for q in eachindex(xi)
            yi[q] = lerp_(x,y,xi[q])
        end
    end 
    return yi
end

function lerp_(x::Union{T,Vector{T}, AbstractRange{T}},y,xi::T) where T <: Real
    j = findfirst(x.>xi)
    
    if isnothing(j) # Deal with extrapolation at the end
        j = length(x)
        i = j-1 
    elseif isone(j) # Deal with extrapolation at the start
        i=1
        j=2
    else
        i = j-1
    end    
    w = abs(x[j]-x[i]) 
    t = (xi-x[i])/w
    yi = (1.0-t)*y[i] + t*y[j]
    return yi
end

"""
    dist(V1,V2)

Computes n-dimensional Euclidean distances

# Description

Function compute an nxm distance matrix for the n inputs points in `V1`, and the
m input points in `V2`. The input points may be multidimensional, in fact they
can be any type supported by the `euclidean` function of `Distances.jl`. 
See also: https://github.com/JuliaStats/Distances.jl
"""
function dist(V1,V2)
    D = Matrix{Float64}(undef,length(V1),length(V2))   
    for i in eachindex(V1)
        for j in eachindex(V2)          
            D[i,j] = euclidean(V1[i],V2[j]) # norm(V1[i]-V2[j])       
        end
    end
    return D
end

function dist(V1::Vector{T},v2::T) where T <: AbstractVector
    D = Matrix{Float64}(undef,length(V1),1)   
    for i in eachindex(V1)        
        D[i,1] = euclidean(V1[i],v2)
    end
    return D
end

function dist(v1::T,V2::Vector{T}) where T <: AbstractVector
    D = Matrix{Float64}(undef,1,length(V2))   
    for j in eachindex(V2)        
        D[1,j] = euclidean(v1,V2[j])
    end
    return D
end

"""
    mindist(V1,V2; getIndex=false, skipSelf = false )

Returns nearest point distances 

# Description

Returns the closest point distance for the input points `V1` with respect to the 
input points `V2`. If the optional parameter `getIndex` is set to `true` (`false` 
by default) then this function also returns the indices of the nearest points 
in `V2` for each point in `V1`. For self-distance evaluation, i.e. if the same 
point set is provided twice, then the optional parameter `skipSelf` can be set 
t0 `true` (default is `false`) if "self distances" (e.g. the nth point to the 
nth point) are to be avoided.  
"""
function mindist(V1,V2; getIndex=false, skipSelf = false )
    D = Vector{Float64}(undef,length(V1))
    d = Vector{Float64}(undef,length(V2))
    if getIndex
        I = Vector{Int64}(undef,length(V1))
    end
    for i in eachindex(V1)
        for j in eachindex(V2)
            if skipSelf && i==j
                d[j] = Inf
            else
                d[j] = euclidean(V1[i],V2[j]) # norm(V1[i]-V2[j]) 
            end       
        end
        if getIndex
            D[i], I[i] = findmin(d)
        else
            D[i] = minimum(d)
        end
    end
    if getIndex
        return D, I
    else
        return D
    end
end


"""
    unique_dict_index(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any

Returns unique values and indices

# Description

Returns the unique entries in `X` as well as the indices for them. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function unique_dict_index(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Nothing}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    for i in eachindex(X)        
        if sort_entries && length(X[1])>1
            x = sort(X[i])
        else
            x = X[i]
        end
        if !haskey(d, x)
            d[x] = nothing
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
        end
    end
    return xUni, indUnique
end


"""
    unique_dict_index_inverse(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any

Returns unique values, indices, and inverse indices

# Description
    
Returns the unique entries in `X` as well as the indices for them and the 
reverse indices to retrieve the original from the unique entries. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function unique_dict_index_inverse(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    indInverse = Vector{Int64}(undef,length(X)) 
    j=0
    for i in eachindex(X)         
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
            indInverse[i] = j  # Store inverse index          
        else
            indInverse[i]=d[x]
        end
    end
    return xUni, indUnique, indInverse
end


"""
    unique_dict_index_count(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any

Returns unique values, indices, and counts

# Description
    
Returns the unique entries in `X` as well as the indices for them and the counts 
in terms of how often they occured. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function unique_dict_index_count(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    c =  Vector{Int64}() 

    j=0
    for i in eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
            push!(c,1)
        else
            c[d[x]] += 1
        end
    end
    return xUni, indUnique, c
end


"""
    unique_dict_index_inverse_count(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any

Returns unique values, indices, inverse indices, and counts

# Description
    
Returns the unique entries in `X` as well as the indices for them and the reverse 
indices to retrieve the original from the unique entries, and also the counts in 
terms of how often they occured. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function unique_dict_index_inverse_count(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    indInverse = Vector{Int64}(undef,length(X)) 
    c =  Vector{Int64}() 

    j=0
    for i in eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
            indInverse[i] = j  # Store inverse index
            push!(c,1)
        else
            indInverse[i]=d[x]
            c[d[x]] += 1
        end
    end
    return xUni, indUnique, indInverse,c
end


"""
    unique_dict_count(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any

Returns unique values and counts

# Description
    
Returns the unique entries in `X` as well as the counts in terms of how often 
they occured. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function unique_dict_count(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    c = Vector{Int64}()
    j = 0
    for i in eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j += 1 # Index in unique array
            d[x] = j # Store count in dict
            push!(xUni, X[i]) # Grow unique set
            push!(c,1) # Grow counting array
        else
            c[d[x]] += 1
        end
    end
    return xUni, c
end


"""
    unique_dict_inverse(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any 

Returns unique values and inverse indices

# Description

Returns the unique entries in `X` as well as the reverse indices to retrieve the 
original from the unique entries. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function unique_dict_inverse(X::Union{Array{T},Tuple{T}}; sort_entries=false) where T <: Any 
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indInverse = Vector{Int64}(undef,length(X)) 

    j=0
    for i in eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            indInverse[i] = j  # Store inverse index
        else
            indInverse[i]=d[x]
        end
    end
    return xUni, indInverse
end 


"""
    unique_dict(X::AbstractVector{T}) where T <: Real    

Returns unique values, indices, and inverse indices. Uses an OrderedDict.

# Description
    
Returns the unique entries in `X` as well as the indices for them and the reverse 
indices to retrieve the original from the unique entries. 
"""
function unique_dict(X::AbstractVector{T}) where T <: Real    
    d = OrderedDict{T ,Int64}() # Use dict to keep track of used values    
    indUnique = Vector{Int64}()
    indReverse = Vector{Int64}(undef,length(X))
    j=0
    for i in eachindex(X)        
        if !haskey(d, X[i])                                          
            j+=1
            d[X[i]] = j # reverse index in dict            
            push!(indUnique, i)                     
            indReverse[i] = j 
        else
            indReverse[i] = d[X[i]]            
        end        
    end
    return collect(keys(d)), indUnique, indReverse
end


"""
    gunique(X; return_unique=true, return_index=false, return_inverse=false, return_counts=false, sort_entries=false)
    
Returns unique values and allows users to choose if they also want: sorting, indices, inverse indices, and counts. 

# Description

Returns the unique entries in `X`. Depending on the optional parameter choices
the indices for the unique entries, the reverse indices to retrieve the original
from the unique entries, as well as counts in terms of how often they occured, 
can be returned. 
The optional parameter `sort_entries` (default is `false`) can be set to `true`
if each entry in X should be sorted, this is helpful to allow the entry [1,2] to 
be seen as the same as [2,1] for instance.  
"""
function gunique(X; return_unique=true, return_index=false, return_inverse=false, return_counts=false, sort_entries=false)
    # Use required unique function 
    if return_unique==true && return_index==false && return_inverse==false && return_counts==false
        # UNIQUE
        if sort_entries && length(X[1])>1
            return unique(sort.(X))
        else
            return unique(X)
        end
    elseif return_unique==true && return_index==true && return_inverse==false && return_counts==false
        # UNIQUE, INDICES
        return unique_dict_index(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==false && return_inverse==false && return_counts==true
        # UNIQUE, COUNTS
        return unique_dict_count(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==false && return_inverse==true && return_counts==false
        # UNIQUE, INVERSE
        return unique_dict_inverse(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==true && return_inverse==true && return_counts==false
        # UNIQUE, INDICES, INVERSE
        return unique_dict_index_inverse(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==true && return_inverse==false && return_counts==true
        # UNIQUE, INDICES, COUNTS
        return unique_dict_index_count(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==true && return_inverse==true && return_counts==true
        # UNIQUE, INDICES, INVERSE, COUNTS
        return unique_dict_index_inverse_count(X; sort_entries=sort_entries)
    end
end


"""
    unique_simplices(F,V=nothing)

Returns unique simplices (such as faces), independant of node order

# Description
    
Returns the unique simplices in F as well as the indices of the unique simplices
and the reverse indices to retrieve the original faces from the unique faces. 
Entries in F are sorted such that the node order does not matter. 
"""
function unique_simplices(F,V=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,F)) 
    else
        n = length(V)
    end
    virtualFaceIndices = sub2ind(n.*ones(Int64,length(F[1])),sort.(F))    
    _, ind1, ind2 = unique_dict(virtualFaceIndices) 

    # Fn, ind1, ind2 = gunique(F; return_unique=true, return_index=true, return_inverse=true, return_counts=false, sort_entries=true)
    # return Fn, ind1, ind2

    return F[ind1], ind1, ind2
end


"""
    ind2sub(siz,ind)

Converts linear indices to subscript indices. 

# Description

Converts the linear indices in `ind`, for a matrix/array with size `siz`, to the 
equivalent subscript indices.  
"""
function ind2sub(siz::Union{Tuple{Vararg{Int64, N}}, Array{Int64, N}},ind::Union{Int64,Tuple{Vararg{Int64, M}}, Array{Int64, M}}) where N where M
    if !isempty(ind) # Not empty so subscript indices will be derived
        numDim = length(siz) # Number of dimensions from length of size
        k = cumprod(siz) # Cumulative product of size
        m = prod(siz) # Number of elements   
        if any(ind.>m) || any(ind.<1)
            throw(BoundsError("Encountered index value out of valid range 1:$m"))
        end
        if isa(ind,Union{Array,Tuple}) # Potentially multiple indices so loop over them
            A = [ind2sub_(ind_i,numDim,k) for ind_i in ind]
        else # This should be a single integer
            A = ind2sub_(ind,numDim,k)      
        end
    else # Empty so return an empty vector
        A = Vector{Int64}[]
    end
    return A
end

# ind2sub helper function to parse just a single linear index and produce a single subscript index set 
function ind2sub_(ind::Int64,numDim::Int64,k::Union{Int64,Array{Int64, N},Tuple{Vararg{Int64, N}}}) where N
    a = Vector{Int64}(undef,numDim) # Initialise a
    for q in numDim:-1:1   # For all dimensions     
        if isone(q) # First 1st dimension
            a[1] = rem(ind-1,k[1]) + 1        
        else       
            p = rem(ind-1,k[q-1]) + 1 # "previous"
            a[q] = (ind - p)./k[q-1] + 1 # Current        
            ind = p # Set indices as "previous"          
        end            
    end    
    return a
end


"""
    sub2ind(siz,A)

Converts subscript indices to linear indices. 

# Description

Converts the subscript indices in `A`, for a matrix/array with size `siz`, to 
the equivalent linear indices.  
"""
function sub2ind(siz::Union{Tuple{Vararg{Int64, N}}, Array{Int64, N}},A::Union{Vector{Vector{Int64}}, Array{NgonFace{M, Int64}, 1}}) where N where M
    numDim = length(siz)
    k = cumprod([siz[i] for i in eachindex(siz)],dims=1)        
    ind = Vector{Int64}(undef,length(A))
    for i in eachindex(A)        
        a = A[i]
        if length(a)==numDim
            if any(a.>siz) || any(a.<1)
                throw(DomainError(A[i],"Indices in A[$i] exceed bounds implied by size provided"))
            end    
        else                
            throw(DimensionMismatch("Incorrect number of indices in  A[$i], size implies number of indices should be $numDim"))
        end            
        ind[i] = sub2ind_(a,numDim,k)
    end                    
    return ind
end

function sub2ind(siz::Union{Tuple{Vararg{Int64, N}}, Vector{Int64}},A::Vector{Int64})  where N  
    numDim = length(siz)  
    if length(A)==numDim
        if any(A.>siz) || any(A.<1)
            throw(DomainError(A,"Indices in A exceed bounds implied by size provided"))
        end 
    else                
        throw(DimensionMismatch("Incorrect number of indices in  A, size implies number of indices should be $numDim"))
    end 
    return sub2ind(siz,[A])[1]
end

function sub2ind_(a::Union{Tuple{Vararg{Int64, N}}, Vector{Int64},NgonFace{M, Int64}},numDim::Int64,k::Union{Int64,Vector{Int64},Tuple{Vararg{Int64, N}}})  where N where M
    if numDim==1 
        ind = a[1]
    else        
        ind = a[1]
        for i=2:numDim
            ind += (a[i].-1).*k[i-1]
        end 
    end
    return ind
end


"""
    meshedges(F::Array{NgonFace{N,T},1}; unique_only=false) where N where T<:Integer   

Returns a mesh's edges.

# Description

This function returns the edges `E` for the input faces defined by `F`. 
The input `F` can either represent a vector of faces or a 
GeometryBasics.Mesh. The convention is such that for a face referring to the 
nodes 1-2-3-4, the edges are 1-2, 2-3, 3-4, 4-1.   
"""
function meshedges(F::Array{NgonFace{N,T},1}; unique_only=false) where N where T<:Integer        
    E = LineFace{Int64}[]    
    for j1 in 1:N # Loop over each node/point for the current simplex           
        if j1<N
            j2 = j1+1
        else
            j2 = 1
        end            
        for f in F # Loop over each simplex        
            push!(E,(f[j1],f[j2]))            
        end 
    end    
    if unique_only # Remove doubles if requested e.g. 1-2 seen as same as 2-1
        E = gunique(E; sort_entries=true);
    end
    return E
end

function meshedges(M::GeometryBasics.Mesh; unique_only=false)   
    return meshedges(faces(M); unique_only=unique_only)
end


"""
    icosahedron(r=1.0)

Creates an icosahedron mesh. 

# Description

Creates a GeometryBasics.Mesh for an icosahedron with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function icosahedron(r=1.0)

    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(ϕ + 2.0)
    t = ϕ*s

    # Create vertices
    V=Vector{Point{3, Float64}}(undef,12)
    V[1 ] = Point{3, Float64}( 0.0,   -s,  -t)
    V[2 ] = Point{3, Float64}( 0.0,   -s,   t)
    V[3 ] = Point{3, Float64}( 0.0,    s,   t)
    V[4 ] = Point{3, Float64}( 0.0,    s,  -t)
    V[5 ] = Point{3, Float64}(  -s,   -t, 0.0)
    V[6 ] = Point{3, Float64}(  -s,    t, 0.0)
    V[7 ] = Point{3, Float64}(   s,    t, 0.0)
    V[8 ] = Point{3, Float64}(   s,   -t, 0.0)
    V[9 ] = Point{3, Float64}(  -t,  0.0,  -s)
    V[10] = Point{3, Float64}(   t,  0.0,  -s)
    V[11] = Point{3, Float64}(   t,  0.0,   s)
    V[12] = Point{3, Float64}(  -t,  0.0,   s)

    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,20)
    F[1 ] = TriangleFace{Int64}(9,4,1)
    F[2 ] = TriangleFace{Int64}(1,5,9)
    F[3 ] = TriangleFace{Int64}(1,8,5)
    F[4 ] = TriangleFace{Int64}(10,8,1)
    F[5 ] = TriangleFace{Int64}(4,10,1)
    F[6 ] = TriangleFace{Int64}(5,2,12)
    F[7 ] = TriangleFace{Int64}(12,2,3)
    F[8 ] = TriangleFace{Int64}(12,3,6)
    F[9 ] = TriangleFace{Int64}(12,6,9)
    F[10] = TriangleFace{Int64}(12,9,5)
    F[11] = TriangleFace{Int64}(10,7,11)
    F[12] = TriangleFace{Int64}(8,10,11)
    F[13] = TriangleFace{Int64}(2,8,11)
    F[14] = TriangleFace{Int64}(3,2,11)
    F[15] = TriangleFace{Int64}(7,3,11)
    F[16] = TriangleFace{Int64}(2,5,8)
    F[17] = TriangleFace{Int64}(10,4,7)
    F[18] = TriangleFace{Int64}(7,6,3)
    F[19] = TriangleFace{Int64}(6,7,4)
    F[20] = TriangleFace{Int64}(6,4,9)
    
    return GeometryBasics.Mesh(V,F)
end


"""
    octahedron(r=1.0)

Creates an octahedron mesh. 
    
# Description

Creates a GeometryBasics.Mesh for an octahedron with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function octahedron(r=1.0)

    s = r/sqrt(2.0)

    # Create vertices    
    V=Vector{Point{3, Float64}}(undef,6)
    V[1 ] = Point{3, Float64}(   -s,    -s,  0.0)
    V[2 ] = Point{3, Float64}(    s,    -s,  0.0)
    V[3 ] = Point{3, Float64}(    s,     s,  0.0)
    V[4 ] = Point{3, Float64}(   -s,     s,  0.0)
    V[5 ] = Point{3, Float64}(  0.0,   0.0,   -r)
    V[6 ] = Point{3, Float64}(  0.0,   0.0,    r)
    
    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,8)
    F[1 ] = TriangleFace{Int64}(5,2,1)
    F[2 ] = TriangleFace{Int64}(5,3,2)
    F[3 ] = TriangleFace{Int64}(5,4,3)
    F[4 ] = TriangleFace{Int64}(5,1,4)
    F[5 ] = TriangleFace{Int64}(6,1,2)
    F[6 ] = TriangleFace{Int64}(6,2,3)
    F[7 ] = TriangleFace{Int64}(6,3,4)
    F[8 ] = TriangleFace{Int64}(6,4,1)
    
    return GeometryBasics.Mesh(V,F)
end


"""
    dodecahedron(r=1.0)

Creates a dodecahedron mesh. 
    
# Description

Creates a GeometryBasics.Mesh for an dodecahedron with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function dodecahedron(r=1.0)

    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(3.0)
    t = ϕ*s    
    w = (ϕ-1.0)*s

    # Create vertices    
    V=Vector{Point{3, Float64}}(undef,20)
    V[1 ] = Point{3, Float64}(   s,   s,   s)
    V[2 ] = Point{3, Float64}(   w, 0.0,   t)
    V[3 ] = Point{3, Float64}(  -t,  -w, 0.0)
    V[4 ] = Point{3, Float64}(   t,   w, 0.0)
    V[5 ] = Point{3, Float64}(  -s,   s,  -s)
    V[6 ] = Point{3, Float64}( 0.0,  -t,  -w)
    V[7 ] = Point{3, Float64}(  -t,   w, 0.0)
    V[8 ] = Point{3, Float64}(   s,  -s,   s)
    V[9 ] = Point{3, Float64}(  -s,   s,   s)
    V[10] = Point{3, Float64}(  -s,  -s,   s)
    V[11] = Point{3, Float64}(   s,  -s,  -s)
    V[12] = Point{3, Float64}(   w, 0.0,  -t)
    V[13] = Point{3, Float64}(  -s,  -s,  -s)
    V[14] = Point{3, Float64}( 0.0,  -t,   w)
    V[15] = Point{3, Float64}( 0.0,   t,  -w)
    V[16] = Point{3, Float64}(  -w, 0.0,   t)
    V[17] = Point{3, Float64}(   t,  -w, 0.0)
    V[18] = Point{3, Float64}(  -w, 0.0,  -t)
    V[19] = Point{3, Float64}(   s,   s,  -s)
    V[20] = Point{3, Float64}( 0.0,   t,   w)

    # Create faces
    F = Vector{NgonFace{5,Int64}}(undef,12)
    F[1 ] = NgonFace{5,Int64}(20, 9,16, 2, 1)
    F[2 ] = NgonFace{5,Int64}( 2,16,10,14, 8)
    F[3 ] = NgonFace{5,Int64}(16, 9, 7, 3,10)
    F[4 ] = NgonFace{5,Int64}( 7, 9,20,15, 5)
    F[5 ] = NgonFace{5,Int64}(18,13, 3, 7, 5)
    F[6 ] = NgonFace{5,Int64}( 3,13, 6,14,10)
    F[7 ] = NgonFace{5,Int64}( 6,13,18,12,11)
    F[8 ] = NgonFace{5,Int64}( 6,11,17, 8,14)
    F[9 ] = NgonFace{5,Int64}(11,12,19, 4,17)
    F[10] = NgonFace{5,Int64}( 1, 2, 8,17, 4)
    F[11] = NgonFace{5,Int64}( 1, 4,19,15,20)
    F[12] = NgonFace{5,Int64}(12,18, 5,15,19)
    
    return GeometryBasics.Mesh(V,F)
end


"""
    cube(r=1.0)

Creates a cube mesh. 

# Description

Creates a GeometryBasics.Mesh for an cube with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function cube(r=1.0)

    # Create vertices       
    s = r/sqrt(3.0)

    V = Vector{Point{3, Float64}}(undef,8)
    V[1 ] = Point{3, Float64}( -s, -s, -s)
    V[2 ] = Point{3, Float64}( -s,  s, -s)
    V[3 ] = Point{3, Float64}(  s,  s, -s)
    V[4 ] = Point{3, Float64}(  s, -s, -s)
    V[5 ] = Point{3, Float64}( -s, -s,  s)
    V[6 ] = Point{3, Float64}( -s,  s,  s)
    V[7 ] = Point{3, Float64}(  s,  s,  s)
    V[8 ] = Point{3, Float64}(  s, -s,  s)
        
    # Create faces
    F = Vector{QuadFace{Int64}}(undef,6)
    F[1 ] = QuadFace{Int64}(1,2,3,4)
    F[2 ] = QuadFace{Int64}(8,7,6,5)
    F[3 ] = QuadFace{Int64}(5,6,2,1)
    F[4 ] = QuadFace{Int64}(6,7,3,2)    
    F[5 ] = QuadFace{Int64}(7,8,4,3)    
    F[6 ] = QuadFace{Int64}(8,5,1,4)    

    return GeometryBasics.Mesh(V,F)
end


"""
    tetrahedron(r=1.0)

Creates a tetrahedron mesh. 

# Description

Creates a GeometryBasics.Mesh for an tetrahedron with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function tetrahedron(r=1.0)

    # Create vertices    
    a = r*sqrt(2.0)/sqrt(3.0)
    b = -r*sqrt(2.0)/3.0
    c = -r/3.0       

    V=Vector{Point{3, Float64}}(undef,4)
    V[1 ] = Point{3, Float64}(   -a,      b,   c)
    V[2 ] = Point{3, Float64}(    a,      b,   c)    
    V[3 ] = Point{3, Float64}(  0.0,    0.0,   r)
    V[4 ] = Point{3, Float64}(  0.0, -2.0*b,   c)  

    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,4)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    F[2 ] = TriangleFace{Int64}(4,2,1)
    F[3 ] = TriangleFace{Int64}(4,3,2)
    F[4 ] = TriangleFace{Int64}(4,1,3)
    
    return GeometryBasics.Mesh(V,F)
end


"""
    platonicsolid(n,r=1.0)

Returns a platonic solid mesh.

# Description

Creates a GeometryBasics mesh description for a platonic solid of choice. The 
input `n` defines the choice.
1. tetrahedron
2. cube
3. octahedron
4. icosahedron
5. dodecahedron

The final input parameter `r` defines the radius of the platonic solid (the 
radius of the circumsphere to the vertices).  The default radius, when not 
supplied, is `1.0`. 

# Arguments

n::Integer, defining platonic solid type
r::Float64, defining circumsphere radius

"""
function platonicsolid(n::Integer,r=1.0)
    if isone(n)
        M = tetrahedron(r)
    elseif n==2
        M = cube(r)
    elseif n==3
        M = octahedron(r)
    elseif n==4 
        M = icosahedron(r)
    elseif n==5
        M = dodecahedron(r)
    end
    return M
end

"""
    tofaces(FM::Vector{Vector{TF}}) where TF<:Integer
    tofaces(FM::Matrix{TF})  where TF<:Integer
    tofaces(FM::Vector{NgonFace{m, OffsetInteger{-1, TF}}} ) where m where TF <: Integer
    tofaces(FM::Vector{NgonFace{m, TF}} ) where m where TF <: Integer 

Converts input to GeometryBasics compliant faces with standard integer types. 

# Description 

The `tofaces` function converts "non-standard" (for Comodo) face set 
descriptions to "standard" ones. The following is considered such as standard: 
`Vector{GeometryBasics.NgonFace{N,T}} where N where T<:Integer` 
The input faces `FM` are converted to this format. `FM` can be of the following 
types: 
* `FM::Vector{Vector{TF}} where TF<:Integer`, whereby each Vector entry is 
considered a face
* `FM::Matrix{TF} where TF<:Integer`, whereby each row is considered a face
* `Vector{NgonFace{m, OffsetInteger{-1, TF}}} where TF<:Integer`, whereby the 
special integer type `OffsetInteger{-1, TF}` is converted to `Int64`.  
If the intput is already of the right type this function leaves the input 
unchanged.
"""
function tofaces(FM::Vector{Vector{TF}}) where TF<:Integer
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)    
    m = length(FM[1]) # Get number of points from first
    if m == 2 # Edges
        F = [LineFace{TF}(f) for f in FM]
    elseif m == 3 # Triangles
        F = [TriangleFace{TF}(f) for f in FM]
    elseif m == 4 # Quads
        F = [QuadFace{TF}(f) for f in FM]
    else # Other mesh type
        F = [NgonFace{m,TF}(f) for f in FM]
    end
    return F
end

function tofaces(FM::Matrix{TF})  where TF<:Integer
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)
    m = size(FM,2) # number of points per face
    if m == 2 # Edges
        F = [LineFace{TF}(f) for f in eachrow(FM)]
    elseif m == 3 # Triangles
        F = [TriangleFace{TF}(f) for f in eachrow(FM)]
    elseif m == 4 # Quads
        F = [QuadFace{TF}(f) for f in eachrow(FM)]        
    else # Other mesh type
        F = [NgonFace{m,TF}(f) for f in eachrow(FM)]
    end
    return F
end

function tofaces(FM::Vector{NgonFace{m, OffsetInteger{-1, TF}}} ) where m where TF <: Integer
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)    
    if m == 2 # Edges
        F = [LineFace{Int64}(f) for f in FM]
    elseif m == 3 # Triangles
        F = [TriangleFace{Int64}(f) for f in FM]
    elseif m == 4 # Quads
        F = [QuadFace{Int64}(f) for f in FM]
    else # Other mesh type
        F = [NgonFace{m,Int64}(f) for f in FM]
    end
    return F
end

function tofaces(FM::Vector{NgonFace{m, TF}} ) where m where TF <: Integer    
    return FM
end

"""
    topoints(VM::Matrix{T}) where T<: Real
    topoints(VM::Union{Array{Vec{N, T}, 1}, GeometryBasics.StructArray{TT,1} }) where TT <: AbstractPoint{N,T} where T <: Real where N   
    topoints(VM::Vector{Vector{T}}) where T <: Real  
    topoints(VM::Vector{Point{ND,TV}}) where ND where TV <: Real    

Converts input to GeometryBasics compliant simple points without meta content.

# Description 
The `topoints` function converts the "non-standard" (for Comodo) input points 
defined by `VM` to the "standard" format:
`VM::Vector{Point{ND,TV}} where ND where TV <: Real`.
For matrix input each row is considered a point. For vector input each vector 
entry is considered a point.     
"""
function topoints(VM::Matrix{T}) where T<: Real
    m = size(VM,2)
    return [Point{m, T}(v) for v in eachrow(VM)]
end

function topoints(VM::Union{Array{Vec{N, T}, 1}, GeometryBasics.StructArray{TT,1} }) where TT <: AbstractPoint{N,T} where T <: Real where N   
    if eltype(VM)<:PointMeta{N,T} where N where T<:Real               
        return VM.position
    else        
        return [Point{N, T}(v) for v in VM]
    end
end

function topoints(VM::Vector{Vector{T}}) where T <: Real    
        m = length(VM[1])
        return [Point{m, T}(v) for v in VM]
end

function topoints(VM::Vector{Point{ND,TV}}) where ND where TV <: Real        
    return VM
end

"""
    togeometrybasics_mesh

Converts the input to a GeometryBasics.Mesh

# Description

This function converts the input faces `F` and vertices `V` to a 
GeometryBasics.Mesh. The function `tofaces` and `topoints` are used prior to 
conversion, to ensure standard faces and point types are used. 
"""
function togeometrybasics_mesh(VM,FM)
    V = topoints(VM)
    F = tofaces(FM)
    return GeometryBasics.Mesh(V,F)
end


"""
    edgecrossproduct(F,V::Vector{Point{ND,T}}) where ND where T<:Real  
    edgecrossproduct(M::GeometryBasics.Mesh)

Returns the edge cross product, useful for nomal direction and area computations. 

    # Description

This function computes the so-called edge-cross-product for a input mesh that is
either defined by the faces `F` and vertices `V` or the mesh `M`. 
"""
function edgecrossproduct(F,V::Vector{Point{ND,TV}}) where ND where TV<:Real
    C = Vector{GeometryBasics.Vec{ND, TV}}(undef,length(F)) # Allocate array cross-product vectors
    n =  length(F[1]) # Number of nodes per face    
    for q in eachindex(F) # Loop over all faces
        c  = cross(V[F[q][n]],V[F[q][1]]) # Initialise as cross product of last and first vertex position vector
        @inbounds for qe in 1:(n-1) # Loop from first to end-1            
            c  += cross(V[F[q][qe]],V[F[q][qe+1]]) # Add next edge contribution          
        end
        C[q] = c./2 # Length = face area, direction is along normal vector
    end
    return C
end

function edgecrossproduct(M::GeometryBasics.Mesh) 
    return edgecrossproduct(faces(M),coordinates(M))
end


"""
    facenormal(F,V; weighting=:area)

Returns the normal directions for each face.

# Description

This function computes the per face normal directions for the input mesh defined 
either by the faces `F` and vertices `V` or by the GeometryBasics mesh `M`. 
"""
function facenormal(F,V::Vector{Point{ND,TV}}) where ND where TV<:Real
    C = edgecrossproduct(F,V)
    return C./norm.(C)
end

function facenormal(M::GeometryBasics.Mesh) 
    return facenormal(faces(M),coordinates(M))
end


"""
    facearea(F,V)

Returns the area for each face. 

# Description

This function computes the per face area for the input mesh defined either by 
the faces `F` and vertices `V` or by the GeometryBasics mesh `M`. 
"""
function facearea(F,V::Vector{Point{ND,TV}}) where ND where TV<:Real
    return norm.(edgecrossproduct(F,V))
end

function facearea(M::GeometryBasics.Mesh)         
    return facearea(faces(M),coordinates(M))
end


"""
    vertexnormal(F,V; weighting=:area)

Returns the surface normal at each vertex.

# Description

This function computes the per vertex surface normal directions for the input 
mesh defined either by the faces `F` and vertices `V` or by the GeometryBasics
mesh `M`. The optional parameter `weighting` sets how the face normal directions 
are averaged onto the vertices. If `weighting=:none` a plain average for the 
surrounding faces is used. If instead `weighting=:area` (default), then the
average is weighted based on the face areas. 
"""
function vertexnormal(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; weighting=:area) where N where TF<:Integer where ND where TV<:Real       
    return normalizevector(simplex2vertexdata(F,facenormal(F,V),V; weighting=weighting))
end

function vertexnormal(M::GeometryBasics.Mesh; weighting=:area)     
    return normalizevector(vertexnormal(faces(M),coordinates(M); weighting=weighting))
end


"""
    edgelengths(E::LineFace,V)
    edgelengths(F,V)
    edgelengths(M::GeometryBasics.Mesh)            

Returns edge lengths.

# Description
This function computes the lengths of the edges defined by edge vector `E` (e.g
as obtained from `meshedges(F,V)`, where `F` is a face vector, and `V` is a 
vector of vertices. 
Alternatively the input mesh can be a GeometryBasics mesh `M`.
"""
function edgelengths(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real
    if eltype(F)<:LineFace{T} where T<:Integer # Already edges 
        return [norm(V[e[1]]-V[e[2]]) for e in F]
    else # Need to compute edges
        return edgelengths(meshedges(F; unique_only=true),V)
    end
end

function edgelengths(M::GeometryBasics.Mesh)        
    return edgelengths(faces(M),coordinates(M))    
end


"""
    subtri(F,V,n; method = :linear)
    subtri(F,V,n; method = :Loop)

Refines triangulations through splitting.

# Description

The `subtri` function refines triangulated meshes iteratively. For each iteration
each original input triangle is split into 4 triangles to form the refined mesh 
(one central one, and 3 at each corner). The following refinement methods are 
implemented: 
    
`method=:linear` : This is the default method, and refines the triangles in a 
simple linear manor through splitting. Each input edge simply obtains a new 
mid-edge node. 
    
`method=:Loop` : This method features Loop-subdivision [1,2]. Rather than linearly 
splitting edges and maintaining the original coordinates, as for the linear 
method, this method computes the new points in a special weighted sense such 
that the surface effectively approaches a "quartic box spline". Hence this 
method both refines and smoothes the geometry through spline approximation. 

# References
1. [Charles Loop, _Smooth Subdivision Surfaces Based on Triangles_, M.S. Mathematics Thesis, University of Utah. 1987.](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/thesis-10.pdf)
2. [Jos Stam, Charles Loop, _Quad/Triangle Subdivision_, doi: 10.1111/1467-8659.t01-2-00647](https://doi.org/10.1111/1467-8659.t01-2-00647)
"""
function subtri(F::Vector{NgonFace{3,TF}},V::Vector{Point{ND,TV}},n::Int64; method = :linear) where TF<:Integer where ND where TV <: Real
    
    if iszero(n)
        return F,V
    elseif isone(n)
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        Fm1 = [TriangleFace{TF}(a.+length(V)) for a in eachrow(reshape(indReverse,length(F),length(F[1])))] 
        Fm2 = Vector{TriangleFace{TF}}(undef,length(Fm1))
        Fm3 = Vector{TriangleFace{TF}}(undef,length(Fm1))
        Fm4 = Vector{TriangleFace{TF}}(undef,length(Fm1))        
        for i in eachindex(F)                        
            Fm2[i] = TriangleFace{TF}([Fm1[i][1], Fm1[i][3], F[i][1]])
            Fm3[i] = TriangleFace{TF}([Fm1[i][2], Fm1[i][1], F[i][2]])
            Fm4[i] = TriangleFace{TF}([Fm1[i][3], Fm1[i][2], F[i][3]])
        end

        # Create combined face set
        Fn = [Fm1; Fm2; Fm3; Fm4]        

        con_E2F = con_edge_face(F,Eu,indReverse)

        # Create new vertices depending on method
        if method == :linear # Simple linear splitting
            # Create complete point set
            Vn = [V; simplexcenter(Eu,V)]  # Old and new mid-edge points          
        elseif method == :Loop #Loop subdivision 
    
            # New mid-edge like vertices
            Vm = Vector{Point{ND,TV}}(undef,length(Eu)) 
            for q in eachindex(Eu) # For each edge index                        
                F_touch = F[con_E2F[q]] # Faces sharing current edge, mostly 2 but 1 for a boundary edge
                indVerticesTouch = Vector{TF}() 
                for f in F_touch        
                    b = f.!=Eu[q][1] .&& f.!=Eu[q][2]      
                    if any(b)  
                        append!(indVerticesTouch,f[b])           
                    end
                end        
                Vm[q]=3/8 .*(V[Eu[q][1]] .+ V[Eu[q][2]])  .+ 1/8 .* (V[indVerticesTouch[1]] .+ V[indVerticesTouch[2]])
            end
    
            # Modified vertices for original vertices
            Vv = Vector{Point{ND,TV}}(undef,length(V))
            for q in eachindex(V)            
                B_vert_face = [any(f.==q) for f in F]
                F_touch = F[B_vert_face] # Faces mostly 2 but 1 for a boundary edge
                indVerticesTouch = Vector{TF}()
                for f in F_touch                
                    indTouch = f[f.!=q]        
                    for i in indTouch 
                        if i ∉ indVerticesTouch 
                            push!(indVerticesTouch,i)
                        end
                    end
                end
                N = length(indVerticesTouch)                
                v_sum = sum(V[indVerticesTouch],dims=1)[1]                
                β = 1/N * (5/8-(3/8 +1/4*cos((2*π)/N))^2)        
                Vv[q] = (1-N*β) .* V[q] .+ β*v_sum                
            end    
            # Create complete point set
            Vn = [Vv;Vm] # Updated orignals and new "mid-edge-ish" points
        else
            throw(ArgumentError("Incorrect method :$method. Use :linear or :loop"))
        end

        return Fn,Vn    

    elseif n>1
        for _ =1:n
            F,V = subtri(F,V,1; method=method)
        end
        return F,V
    else
        throw(ArgumentError("n should be larger than 0"))
    end
end

"""
subquad(F::Vector{NgonFace{4,TF}},V::Vector{Point{ND,TV}},n::Int64; method=:linear) where TF<:Integer where ND where TV <: Real
subquad(F::Vector{NgonFace{4,TF}},V::Vector{Point{ND,TV}},n::Int64; method=:Catmull_Clark) where TF<:Integer where ND where TV <: Real

Refines quadrangulations through splitting.

# Description

The `subquad` function refines quad meshes iteratively. For each iteration each
original input quad is split into 4 smaller quads to form the refined mesh. 
The following refinement methods are implemented: 
    
`method=:linear` : This is the default method, and refines the quads in a 
simple linear manor through splitting. Each input edge simply obtains a new 
mid-edge node, and each face obtains a new central node. 
    
`method=:Catmull_Clark` : This method features Catmull_Clark-subdivision [1]. 
Rather than linearly splitting edges and maintaining the original coordinates, 
as for the linear method, this method computes the new points in a special 
weighted sense such that the surface effectively approaches a bicubic B-spline 
surface. Hence this method both refines and smoothes the geometry through 
spline approximation. 

# References
1. [E. Catmull and J. Clark, _Recursively generated B-spline surfaces on arbitrary topological meshes_, Computer-Aided Design, vol. 10, no. 6, pp. 350-355, Nov. 1978, doi: 10.1016/0010-4485(78)90110-0](https://doi.org/10.1016/0010-4485(78)90110-0).
"""
function subquad(F::Vector{NgonFace{4,TF}},V::Vector{Point{ND,TV}},n::Int64; method=:linear) where TF<:Integer where ND where TV <: Real
    if iszero(n)
        return F,V
    elseif isone(n)

        # Get edges
        E = meshedges(F) # Non-unique edges

        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        con_F2E = con_face_edge(F,Eu,indReverse)
        con_E2F = con_edge_face(F,Eu,indReverse)
        con_V2E = con_vertex_edge(Eu,V)
        con_V2F = con_vertex_face(F,V)

        # Define vertices
        if method ==:linear
            Ve = simplexcenter(Eu,V) # Mid edge points
            Vf = simplexcenter(F,V)  # Mid face points
            Vn = [V;Ve;Vf] # Joined point set
        elseif method ==:Catmull_Clark
            # Mid face points
            Vf = simplexcenter(F,V)  
            Ve_mid = simplexcenter(Eu,V) # Mid edge points

            # Edge points 
            Ve = Vector{Point{ND,TV}}(undef,length(Eu)) # Initialize edge points
            for q in eachindex(Eu)                         
                Ve[q] = (mean(Vf[con_E2F[q]],dims=1)[1] .+ Ve_mid[q])./2.0
            end

            # Vertex points 
            Vv = Vector{Point{ND,TV}}(undef,length(V)) # Initialize vertex points
            for q in eachindex(V) # Loop over all vertices
                indF = con_V2F[q]
                indE = con_V2E[q]
                N = length(indF) # Number of faces (or edges) touching this vertex                    
                Vv[q] = (mean(Vf[indF],dims=1)[1] .+ 2.0.*mean(Ve_mid[indE],dims=1)[1] .+ (N-3.0).*V[q])./N
            end

            Vn = [Vv;Ve;Vf] # Joined point set
        else
            throw(ArgumentError("Incorrect method :$method. Use :linear or :Catmull_Clark"))
        end

        # Define faces
        Fn = Vector{QuadFace{TF}}(undef,length(F)*4)        
        nv = length(V)
        ne = length(Eu)
        for q in eachindex(F)
            i = 1 + (q-1)*4
            for ii = 0:3
                Fn[i+ii] = QuadFace{TF}([F[q][ii+1], con_F2E[q][ii+1]+nv, q+nv+ne, con_F2E[q][1+mod(3+ii,4)]+nv])                
            end            
        end
        return Fn,Vn
    elseif n>1
        for _ =1:n
            F,V = subquad(F,V,1;method=method)
        end
        return F,V
    else
        throw(ArgumentError("n should be larger than 0"))
    end
end


"""
    geosphere(n::Int64,r::T) where T <: Real

Returns a geodesic sphere triangulation

# Description

This function returns a geodesic sphere triangulation based on the number of
refinement iterations `n` and the radius `r`. Geodesic spheres (aka Buckminster-Fuller
 spheres) are triangulations of a sphere that have near uniform edge lenghts. 
The algorithm starts with a regular icosahedron. Next this icosahedron is refined 
`n` times, while nodes are pushed to a sphere surface with radius `r` at each
iteration. 
"""
function geosphere(n::Int64,r::T) where T <: Real
    M = platonicsolid(4,r)
    V = coordinates(M)
    F = faces(M)
    for _ = 1:n
        F,V = subtri(F,V,1)
        for q in eachindex(V)
            v = V[q]
            rn = sqrt(sum(v.^2))
            V[q] = v .* (r/rn)
        end
    end
    return F,V
end

"""
    hexbox(boxDim,boxEl)

Returns a hexahedral mesh of a box

# Description

This function returns a hexahedral mesh for a 3D rectangular box domain. 
"""
function hexbox(boxDim,boxEl)
    boxNod = boxEl.+1 # Number of nodes in each direction
    numElements = prod(boxEl) # Total number of elements
    numNodes = prod(boxNod) # Total number of nodes

    # Create hexahedral element description 
    ind1 = 1:numElements
    ijk_shift = [[0,0,0], 
                [1,0,0],
                [1,1,0],
                [0,1,0],
                [0,0,1],
                [1,0,1],
                [1,1,1],
                [0,1,1]]
                
    E = [Vector{Int64}(undef,8) for _ in 1:numElements] # Allocate elements

    @inbounds for q in 1:numElements
        ijk_1 = ind2sub(boxEl,ind1[q])    
        ijk_2 = ijk_1 .+ ijk_shift[2]
        ijk_3 = ijk_1 .+ ijk_shift[3]
        ijk_4 = ijk_1 .+ ijk_shift[4]
        ijk_5 = ijk_1 .+ ijk_shift[5]
        ijk_6 = ijk_1 .+ ijk_shift[6]
        ijk_7 = ijk_1 .+ ijk_shift[7]
        ijk_8 = ijk_1 .+ ijk_shift[8]
        
        E[q] = [sub2ind(boxNod,ijk_1),sub2ind(boxNod,ijk_2),sub2ind(boxNod,ijk_3),sub2ind(boxNod,ijk_4),sub2ind(boxNod,ijk_5),sub2ind(boxNod,ijk_6),sub2ind(boxNod,ijk_7),sub2ind(boxNod,ijk_8)]
    end

    # Create vertices aka nodal coordinates
    indNodes = collect(Int64,1:numNodes)
    IJK_nodes = ind2sub(boxNod,indNodes)

    V = convert(Vector{Point{3, Float64}},IJK_nodes)
    for q in eachindex(V)
        V[q]=(V[q].-[1,1,1]).*(boxDim./boxEl)
    end

    # Create face sets from elements
    F = Vector{QuadFace{Int64}}(undef,numElements*6) # Allocate faces
    CF_type = Vector{Int64}(undef,numElements*6) # Allocate face color/label data
    for q = eachindex(E) # Loop over all elements
        e = E[q] # The current element 
        qf = 1 + (q-1)*6 # Index mapping into face array 

        # Add the current element's faces
        F[qf  ] = QuadFace{Int64}(e[4],e[3],e[2],e[1]) # Top
        F[qf+1] = QuadFace{Int64}(e[5],e[6],e[7],e[8]) # Bottom
        F[qf+2] = QuadFace{Int64}(e[1],e[2],e[6],e[5]) # Side 1
        F[qf+3] = QuadFace{Int64}(e[3],e[4],e[8],e[7]) # Side 2
        F[qf+4] = QuadFace{Int64}(e[2],e[3],e[7],e[6]) # Front
        F[qf+5] = QuadFace{Int64}(e[4],e[1],e[5],e[8]) # Back

        # Add the current element's face color/labels
        CF_type[qf  ] = 1 # Top
        CF_type[qf+1] = 2 # Bottom
        CF_type[qf+2] = 3 # Side 1
        CF_type[qf+3] = 4 # Side 2
        CF_type[qf+4] = 5 # Front
        CF_type[qf+5] = 6 # Back    
    end

    F_uni,indUni,c_uni = gunique(F,return_index=true, return_inverse=false,return_counts=true,sort_entries=true)
    Lb = isone.(c_uni)
    Fb = F_uni[Lb]
    CF_type_uni = CF_type[indUni]
    CFb_type = CF_type_uni[Lb]

    return E,V,F,Fb,CFb_type
end


function con_face_edge(F,E_uni=nothing,indReverse=nothing)
    if isnothing(E_uni) | isnothing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    return [Vector{Int64}(a) for a in eachrow(reshape(indReverse,length(F),length(F[1])))] # [indReverse[[1,2,3].+ (i-1)*3] for i in eachindex(F)]
end


function con_edge_face(F,E_uni=nothing,indReverse=nothing)
    if isnothing(E_uni) || isnothing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    con_F2E = con_face_edge(F,E_uni,indReverse)
    
    con_E2F = [Vector{Int64}() for _ in 1:length(E_uni)]
    for i_f in eachindex(F)
        for i in con_F2E[i_f]
            push!(con_E2F[i],i_f)
        end
    end 
    return con_E2F
end


function con_face_face(F,E_uni=nothing,indReverse=nothing,con_E2F=nothing,con_F2E=nothing)
    if length(F)>1 # More than one face so compute connectivity
        if isnothing(E_uni)| isnothing(indReverse)
            E = meshedges(F)
            E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        end    
        if isnothing(con_E2F) 
            con_E2F = con_edge_face(F,E_uni)
        end
        if isnothing(con_F2E)
            con_F2E = con_face_edge(F,E_uni,indReverse)                 
        end
        con_F2F = [Vector{Int64}() for _ in 1:length(F)]
        for i_f in eachindex(F)
            for i in reduce(vcat,con_E2F[con_F2E[i_f]])    
                if i!=i_f     
                    push!(con_F2F[i_f],i)
                end 
            end
        end
        return con_F2F
    else # Just one face, so return empty
        return [Vector{Int64}() ]
    end
end


function con_face_face_v(F,V=nothing,con_V2F=nothing)
    if length(F)>1 # More than one face so compute connectivity
        if isnothing(con_V2F) 
            con_V2F = con_vertex_face(F,V)  # VERTEX-FACE connectivity
        end
        con_F2F = [Vector{Int64}() for _ in 1:length(F)]
        for i_f in eachindex(F)
            for i in unique(reduce(vcat,con_V2F[F[i_f]]))    
                if i!=i_f     
                    push!(con_F2F[i_f],i)
                end 
            end
        end
        return con_F2F
    else # Just one face, so return empty
        return [Vector{Int64}() ]
    end
end


function con_vertex_simplex(F,V=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,F))
    else
        n = length(V)
    end
    con_V2F = [Vector{Int64}() for _ in 1:n]
    for i_f in eachindex(F)
        for i in F[i_f]
            push!(con_V2F[i],i_f)
        end
    end
    return con_V2F
end


function con_vertex_face(F,V=nothing)
    return con_vertex_simplex(F,V)
end


function con_vertex_edge(E,V=nothing)
    return con_vertex_simplex(E,V)
end


function con_edge_edge(E_uni,con_V2E=nothing)
    if isnothing(con_V2E)
        con_V2E = con_vertex_edge(E_uni) 
    end    
    con_E2E = [Vector{Int64}() for _ in 1:length(E_uni)]
    for i_e in eachindex(E_uni)
        for i in reduce(vcat,con_V2E[E_uni[i_e]])    
            if i!=i_e     
                push!(con_E2E[i_e],i)
            end 
        end
    end
    return con_E2E
end


function con_vertex_vertex_f(F,V=nothing,con_V2F=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,F))
    else 
        n = length(V)
    end
    
    if isnothing(con_V2F)
        con_V2F = con_vertex_face(F,V)
    end

    con_V2V = [Vector{Int64}() for _ in 1:n]
    @inbounds for i_v in 1:n
        if !isempty(con_V2F[i_v])
            for i in unique(reduce(vcat,F[con_V2F[i_v]]))
                if i_v!=i
                    push!(con_V2V[i_v],i)
                end
            end
        end
    end
    return con_V2V

end


function con_vertex_vertex(E,V=nothing,con_V2E=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,E))
    else 
        n = length(V)
    end
    if isnothing(con_V2E)
        con_V2E = con_vertex_edge(E,V)
    end

    con_V2V = [Vector{Int64}() for _ in 1:n]
    @inbounds for i_v in 1:n
        if !isempty(con_V2E[i_v])
            for i in unique(reduce(vcat,E[con_V2E[i_v]]))
                if i_v!=i
                    push!(con_V2V[i_v],i)
                end
            end
        end
    end

    return con_V2V
end

function meshconnectivity(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real

    # EDGE-VERTEX connectivity
    E = meshedges(F)
    E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

    # FACE-EDGE connectivity
    con_F2E = con_face_edge(F,E_uni,indReverse)    

    # EDGE-FACE connectivity
    con_E2F = con_edge_face(F,E_uni)

    # FACE-FACE connectivity wrt edges
    con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)

    # VERTEX-FACE connectivity
    con_V2F = con_vertex_face(F,V)

    # VERTEX-EDGE connectivity
    con_V2E = con_vertex_edge(E_uni,V)

    # EDGE-EDGE connectivity
    con_E2E = con_edge_edge(E_uni,con_V2E)

    # VERTEX-VERTEX connectivity wrt edges
    con_V2V = con_vertex_vertex(E_uni,V,con_V2E)

    # VERTEX-VERTEX connectivity wrt faces
    con_V2V_f = con_vertex_vertex_f(E_uni,V,con_V2E)

    # FACE-FACE connectivity wrt vertices
    con_F2F_v = con_face_face_v(F,con_V2F)

    return ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V, con_V2V_f, con_F2F_v) 
end 

function mergevertices(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; roundVertices = true, numDigitsMerge=nothing) where N where TF<:Integer where ND where TV<:Real

    m = length(V)
    if roundVertices
        if isnothing(numDigitsMerge)
            E = meshedges(F)
            d = [sqrt( sum((V[e[1]] .- V[e[2]]).^2) ) for e in E]
            pointSpacing = mean(d)            
            numDigitsMerge = 6-round(Int64,log10(pointSpacing))
        end

        # Create rounded coordinates to help obtain unique set
        VR = [round.(v,digits = numDigitsMerge) for v in V]

        # Get unique indices and reverse for rounded vertices
        _,ind1,ind2 = gunique(VR; return_index=true, return_inverse=true,sort_entries=false)
        V = V[ind1] # The unique node set
    else
        V,ind1,ind2 = gunique(V; return_index=true, return_inverse=true,sort_entries=false)
    end
    
    if length(V) != m # If the length has changed
        # Correct indices for faces
        for q in eachindex(F)
            F[q] = ind2[F[q]]
        end
    end

    return F,V,ind1,ind2
end

"""

    smoothmesh_laplacian(F,V,con_V2V=nothing; n=1, λ=0.5)

# Description

This function implements weighted Laplacian mesh smoothing. At each 
iteration, this method replaces each point by an updated coordinate based on the 
mean coordinates of that point's Laplacian umbrella. The update features a lerp
like weighting between the previous iterations coordinates and the mean 
coordinates. The code features `Vs[q] = (1.0-λ).*Vs[q] .+ λ*mean(V[con_V2V[q]])`
As can be seen, the weighting is controlled by the input parameter `λ` which is
in the range (0,1). If `λ=0` then no smoothing occurs. If `λ=1` then pure 
Laplacian mean based smoothing occurs. For intermediate values a linear blending
between the two occurs.  
"""
function smoothmesh_laplacian(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}, n=1, λ=0.5; con_V2V=nothing, constrained_points=nothing) where N where TF<:Integer where ND where TV<:Real
    
    if λ>1.0 || λ<0.0
        throw(ArgumentError("λ should be in the range 0-1"))
    end
    
    if λ>0.0
        if n==0
            return V
        elseif n>0
            # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if nothing
            if isnothing(con_V2V)
                E_uni = meshedges(F;unique_only=true)
                con_V2V = con_vertex_vertex(E_uni)
            end        
            for _ = 1:n
                Vs = deepcopy(V)
                for q in eachindex(V)                
                    Vs[q] = (1.0-λ).*Vs[q] .+ λ*mean(V[con_V2V[q]]) # Linear blend between original and pure Laplacian
                end
                if !isnothing(constrained_points)
                    Vs[constrained_points] = V[constrained_points] # Put back constrained points
                end
                V = Vs
            end
        else #n<0
            throw(ArgumentError("n should be greater or equal to 0"))
        end
    end
    return V
end


"""
    smoothmesh_hc(F,V, con_V2V=nothing; n=1, α=0.1, β=0.5, tolDist=nothing)

# Description 

This function implements HC (Humphrey's Classes) smoothing [1]. This method uses
Laplacian like smoothing but aims to compensate for shrinkage/swelling by also 
"pushing back" towards the original coordinates. 

# Reference 
1. [Vollmer et al., _Improved Laplacian Smoothing of Noisy Surface Meshes_, 1999. doi: 10.1111/1467-8659.00334](https://doi.org/10.1111/1467-8659.00334)
"""
function smoothmesh_hc(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}, n=1, α=0.1, β=0.5; con_V2V=nothing, tolDist=nothing, constrained_points=nothing) where N where TF<:Integer where ND where TV<:Real

    if α>1.0 || α<0.0
        throw(ArgumentError("α should be in the range 0-1"))
    end

    if β>1.0 || β<0.0
        throw(ArgumentError("β should be in the range 0-1"))
    end

    if n<0
        throw(ArgumentError("n should greater or equal to 0"))
    end

    if n == 0 
        return V
    else
        # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if nothing
        if isnothing(con_V2V)
            E = meshedges(F)
            E_uni,_ = gunique(E; return_unique=true, return_index=true, return_inverse=false, sort_entries=true)    
            # E_uni,_,_ = unique_simplices(E)
            con_V2V = con_vertex_vertex(E_uni)
        end
        P = deepcopy(V) # Copy original input points
        B = deepcopy(V) # Initialise B
        c = 0
        while c<n       
            Q = deepcopy(P) # Reset Q as P for this iteration
            for i in eachindex(V)
                P[i] = mean(Q[con_V2V[i]]) # Laplacian 
                # Compute different vector between P and a point between original 
                # point and Q (which is P before laplacian)
                B[i] = P[i] .- (α.*V[i] .+ (1.0-α).*Q[i])
            end
            d = 0.0        
            for i in eachindex(V)      
                # Push points back based on blending between pure difference vector
                # B and the Laplacian mean of these      
                P[i] = P[i] .- (β.*B[i] .+ (1.0-β).* mean(B[con_V2V[i]]))
            end   
            c+=1 
            if !isnothing(tolDist) # Include tolerance based termination
                d = 0.0
                for i in eachindex(V)
                    d+=sqrt(sum((P[i].-Q[i]).^2)) # Sum of distances
                end
                if d<tolDist # Sum of distance smaller than tolerance?
                    break
                end            
            end
            if !isnothing(constrained_points)
                P[constrained_points] = V[constrained_points] # Put back constrained points
            end
        end
        return P
    end
end

function quadplate(plateDim,plateElem)
    num_x = plateElem[1]+1
    num_y = plateElem[2]+1
    V = gridpoints(range(-plateDim[1]/2,plateDim[1]/2,num_x),range(-plateDim[2]/2,plateDim[2]/2,num_y),0.0)
    # V = Vector{Point{3, Float64}}()
    # for y = range(-plateDim[2]/2,plateDim[2]/2,num_y)
    #     for x = range(-plateDim[1]/2,plateDim[1]/2,num_x)
    #         push!(V,Point{3, Float64}(x,y,0.0))
    #     end
    # end

    F = Vector{QuadFace{Int64}}(undef,prod(plateElem))
    ij2ind(i,j) = i + ((j-1)*num_x) # function to convert subscript to linear indices    
    c = 1
    for i = 1:plateElem[1]
        for j = 1:plateElem[2]        
            F[c] = QuadFace{Int64}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)])
            c += 1
        end
    end
    return F, V
end

function quadsphere(n::Int64,r::T) where T <: Real
    M = platonicsolid(2,r)
    F = faces(M)
    V = coordinates(M)
    if n > 0
        for _ in 1:n
            F,V = subquad(F,V,1)
            V = r .* (V ./ norm.(V))
        end
    end
    return F,V
end

function simplex2vertexdata(F,DF,V=nothing; con_V2F=nothing, weighting=:none)
    
    if isnothing(con_V2F)
        con_V2F = con_vertex_face(F,V)
    end
    if weighting==:area
        if isnothing(V)
           throw(ArgumentError("Vertices need to be provided for area based weighting."))
        else
            A = facearea(F,V)
        end
    end    
    DV = (typeof(DF))(undef,length(con_V2F))
    T = eltype(DV)
    for q in eachindex(DV)
        if weighting==:none
            DV[q] = mean(T,DF[con_V2F[q]])
        elseif weighting==:area            
            a = A[con_V2F[q]]
            DV[q] = sum(T,DF[con_V2F[q]].*a)./sum(a)
        end
    end
    return DV
end

function vertex2simplexdata(F,DV)
    T = eltype(DV) # Element type of data in DV
    DF =  (typeof(DV))(undef,length(F)) # Allocate data for F
    for q in eachindex(F)
        DF[q] = mean(T,DV[F[q]]) # The mean of the vertex data for each entry in F
    end
    return DF
end

function simplexcenter(F,V::Vector{Point{ND,TV}}) where ND where TV<:Real
    return vertex2simplexdata(F,V)
end

function normalizevector(A::Union{Vector{Point{ND,TV}},Vector{Vec{ND,TV}}}) where ND where TV<:Real
        return A./norm.(A)
end

function normalizevector(A::Union{Point{ND,TV},Vec{ND,TV}}) where ND where TV<:Real    
    return A./norm(A)    
end

function circlepoints(r::T,n::Int64; dir=:acw) where T <: Real
    if dir==:acw
        return [Point{3, Float64}(r*cos(t),r*sin(t),0) for t in range(0.0,2.0*π-(2.0*π)/n,n)]
    elseif dir==:cw
        return [Point{3, Float64}(r*cos(t),r*sin(t),0) for t in range(0.0,(2.0*π)/n-2.0*π,n)]
    else
        throw(ArgumentError("Invalid dir specified :$dir, use :acw or :cw"))
    end
end

function circlepoints(f::FunctionType,n::Int64; dir=:acw) where {FunctionType <: Function}
    if dir==:acw
        return [Point{3, Float64}(f(t)*cos(t),f(t)*sin(t),0) for t in range(0,2*π-(2*π)/n,n)]
    elseif dir==:cw
        return [Point{3, Float64}(f(t)*cos(t),f(t)*sin(t),0) for t in range(0,(2*π)/n-2*π,n)]
    end
end

"""
    loftlinear(V1,V2;num_steps=2,close_loop=true,face_type=:tri)

Loft a surface mesh between two input curves

# Description 

The `loftlinear` function spans a surface from input curve `V1` to curve `V2`. 
The surface is formed by "lerping" curves from `V1` to `V2` in `num_loft` 
steps, and forming mesh faces between each curve. If `close_loop==true`
then it is assumed the curves (and therefore the output surface mesh should be 
closed over, i.e. that a connection should be made between each curve end and 
start point. The user can request different face types for the output. The 
default is `face_type=:tri` which will form isoceles triangles (or equilateral 
triangles if the spacing is even) for a planar curve. The other `face_type`
options supported are `:quad` (quadrilateral), and `:tri_slash`. For the 
latter, triangles are formed by slashing the quads.  

# Arguments:
- `V1::Vector`: n-vector 
- `V2::Vector`: n-vector

"""
function loftlinear(V1::Vector{Point{ND,TV}},V2::Vector{Point{ND,TV}};num_steps=nothing,close_loop=true,face_type=:quad) where ND where TV<:Real

    # Derive num_steps from distance and mean curve point spacing if missing    
    if isnothing(num_steps)
        d = mean([norm(V1[i]-V2[i]) for i in eachindex(V1)])
        dp = 0.5* (pointspacingmean(V1)+pointspacingmean(V2))
        num_steps = ceil(Int64,d/dp)        
    end

    # Linearly blending points from first to last
    V = Vector{eltype(V1)}()
    for q in range(0,num_steps,num_steps)
        λ = q/num_steps
        Vn = (1.0-λ).*V1 .+ λ.* V2  
        append!(V,Vn)
    end   
    return loftpoints2surf(V,num_steps;close_loop=close_loop,face_type=face_type) # Return faces and vertices
end 


function loftpoints2surf(V::Vector{Point{ND,TV}},num_steps; close_loop=true,face_type=:quad) where ND where TV<:Real

    # Get number of points in each offset curve
    nc = length(V)/num_steps # Number of points in curve
    if !isinteger(nc) || nc<1
        throw(ArgumentError("The length(V)/num_steps should produce an integer >1 but is instead $nc"))
    end
    nc = Int64(nc)

    # Form faces
    if face_type == :tri
        V0 = deepcopy(V)
        for qq in 2:2:num_steps-1
            i = (1:nc) .+ (qq-1) *nc
            @inbounds for q in 1:nc    
                if q == 1 
                    if close_loop == true      
                        V[i[q]] = 0.5 .* (V0[i[q]]+V0[i[q+1]])
                    end
                elseif q == nc 
                    if close_loop == true      
                        V[i[q]] = 0.5 .* (V0[i[q]]+V0[i[1]]) 
                    end
                else
                    V[i[q]] = 0.5 .* (V0[i[q]]+V0[i[q+1]])
                end
            end
        end 
    end

    ij2ind(i,j) = i + nc*(j-1) # function to convert subscript to linear indices

    # Build faces
    if face_type == :quad    
        F = Vector{QuadFace{Int64}}()
        @inbounds for i = 1:(nc-1)
            @inbounds for j = 1:(num_steps-1)    
                push!(F,QuadFace{Int64}([ij2ind(i,j+1),ij2ind(i+1,j+1),ij2ind(i+1,j),ij2ind(i,j)  ]))
            end
        end

        # Add faces to close over shape if requested
        if close_loop
            @inbounds for q in 1:(num_steps-1)                
                push!(F,QuadFace{Int64}([ ij2ind(nc,q+1), ij2ind(1,q+1), ij2ind(1,q), ij2ind(nc,q) ])) 
            end
        end
    elseif face_type == :tri_slash 
        F = Vector{TriangleFace{Int64}}()
        @inbounds for i = 1:nc-1
            @inbounds for j = 1:num_steps-1    
                push!(F,TriangleFace{Int64}([ ij2ind(i+1,j+1), ij2ind(i+1,j), ij2ind(i,j)     ])) # 1 2 3
                push!(F,TriangleFace{Int64}([ ij2ind(i,j),     ij2ind(i,j+1), ij2ind(i+1,j+1) ])) # 3 4 1
            end
        end

        # Add faces to close over shape if requested
        if close_loop
            @inbounds for q in 1:num_steps-1
                push!(F,TriangleFace{Int64}([ ij2ind(1,q+1), ij2ind(1,q),    ij2ind(nc,q)  ])) # 1 2 3
                push!(F,TriangleFace{Int64}([ ij2ind(nc,q),  ij2ind(nc,q+1), ij2ind(1,q+1) ])) # 3 4 1
            end
        end
    elseif face_type == :tri 
        F = Vector{TriangleFace{Int64}}()
        @inbounds for i = 1:nc-1
            @inbounds for j = 1:(num_steps-1)    
                if iseven(j) # Normal slash
                    push!(F,TriangleFace{Int64}([ ij2ind(i+1,j+1), ij2ind(i+1,j),  ij2ind(i,j)     ])) # 1 2 3
                    push!(F,TriangleFace{Int64}([ ij2ind(i,j),     ij2ind(i,j+1),  ij2ind(i+1,j+1) ])) # 3 4 1
                else # Other slash 
                    push!(F,TriangleFace{Int64}([ ij2ind(i,j+1), ij2ind(i+1,j+1), ij2ind(i+1,j) ])) # 2 3 4
                    push!(F,TriangleFace{Int64}([ ij2ind(i+1,j), ij2ind(i,j),     ij2ind(i,j+1) ])) # 4 1 2                 
                end
            end
        end

        # Add faces to close over shape if requested
        if close_loop
            @inbounds for q in 1:(num_steps-1)
                if iseven(q) 
                    push!(F,TriangleFace{Int64}([ ij2ind(nc,q),  ij2ind(nc,q+1), ij2ind(1,q+1) ])) 
                    push!(F,TriangleFace{Int64}([ ij2ind(1,q+1), ij2ind(1,q),    ij2ind(nc,q)  ])) 
                else
                    push!(F,TriangleFace{Int64}([ ij2ind(nc,q+1), ij2ind(1,q+1), ij2ind(1,q)    ]))
                    push!(F,TriangleFace{Int64}([ ij2ind(1,q),    ij2ind(nc,q),  ij2ind(nc,q+1) ]))   
                end
            end
        end
    elseif face_type ==:quad2tri
        F,V = loftpoints2surf(V,num_steps; close_loop=close_loop,face_type=:quad)
        F = quad2tri(F,V; convert_method = :angle)
    else
        throw(ArgumentError("Invalid face_type specified :$face_type, use :quad, :tri, :tri_slash, or :quad2tri"))
    end
    return F, V
end


function dirplot(ax,V::Vector{Point{ND,TV1}},U::Union{Vector{Point{ND,TV2}},Vector{Vec{ND,TV2}}}; color=:black,linewidth=3,scaleval=1.0,style=:from) where ND where TV1 <: Real where TV2 <: Real
    E = [LineFace{Int64}(i,i+length(V)) for i in 1:length(V)]

    if style==:from        
        P = deepcopy(V)
        append!(P,V.+(scaleval.*U))
    elseif style==:to
        P = V.-(scaleval.*U)
        append!(P,V)        
    elseif style==:through
        UU = (scaleval.*U)/2
        P = V.-UU
        append!(P,V.+UU)        
    else
        throw(ArgumentError("Invalid style specified :$style, use :from, :to, or :through"))
    end    
    hp = wireframe!(ax,GeometryBasics.Mesh(P,E),linewidth=linewidth, transparency=false, color=color)
    return hp
end

function normalplot(ax,F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; type_flag=:face, color=:black,linewidth=3,scaleval=nothing) where N where TF<:Integer where ND where TV<:Real
    if isnothing(scaleval)
        scaleval = pointspacingmean(F,V)/2.0
    end
    if type_flag == :face        
        NF = facenormal(F,V)
        V = simplexcenter(F,V)        
    elseif type_flag == :vertex
        NF = vertexnormal(F,V)          
    else
        throw(ArgumentError("Incorrect type_flag, use :face or :vertex"))
    end 
    hp = dirplot(ax,V,NF; color=color,linewidth=linewidth,scaleval=scaleval,style=:from)
    return hp 
end

function normalplot(ax,M::GeometryBasics.Mesh; type_flag=:face, color=:black,linewidth=3,scaleval=nothing)
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    return normalplot(ax,F,V;type_flag=type_flag, color=color,linewidth=linewidth,scaleval=scaleval)
end

function wrapindex(i::Union{Array{Int64,1},UnitRange{Int64},StepRange{Int64, Int64}}, n)
    return 1 .+ mod.(i .+ (n-1),n)
end 

function wrapindex(i::Int64,n)
    return 1+mod(i+(n-1),n)
end

function edgeangles(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real
    m = length(F[1])
    A = Vector{GeometryBasics.Vec{m, Float64}}()
    for f in F
        a = Vector{Float64}(undef,m)
        for i in 1:m                        
            ip1 = wrapindex(i+1,m)            
            ip2 = wrapindex(i+2,m)
            n1 = normalizevector(V[f[ip1]]-V[f[i]])
            n2 = normalizevector(V[f[ip2]]-V[f[ip1]])
            a[i] = acos(clamp(dot(n1,n2),-1.0,1.0))
        end
        push!(A,a)
    end
    return A
end

function quad2tri(F::Vector{QuadFace{TF}},V::Vector{Point{ND,TV}}; convert_method = :angle) where TF<:Integer where ND where TV<:Real
    # Local functions for slash based conversion   
    forw_slash(f) = [TriangleFace{TF}(f[1],f[2],f[3]),TriangleFace{TF}(f[3],f[4],f[1])] # Forward slash 
    back_slash(f) = [TriangleFace{TF}(f[1],f[2],f[4]),TriangleFace{TF}(f[2],f[3],f[4])] # Back slash

    Ft = Vector{TriangleFace{TF}}()#(undef,length(Fn1)*2)
    for f in F        
        if convert_method == :forward
            ft = forw_slash(f)
        elseif convert_method == :backward
            ft = back_slash(f)
        elseif convert_method == :angle 
            ff = forw_slash(f)
            fb = back_slash(f)
            # Get edge angles
            af = edgeangles(ff,V)
            ab = edgeangles(fb,V)
            # Compare angles to perfect triangle angle π/3
            δaf = maximum(abs.(reduce(vcat,af) .- (π/3)))
            δab = maximum(abs.(reduce(vcat,ab) .- (π/3)))
            if δaf<δab
                ft = ff
            else
                ft = fb
            end    
        else
            throw(ArgumentError("Incorrect convert_method set $convert_method, use :forward, :backward, or :angle"))
        end
        push!(Ft,ft[1])
        push!(Ft,ft[2])
    end
    return Ft
end

function remove_unused_vertices(F,V::Vector{Point{ND,TV}})::Tuple where ND where TV<:Real
    if isempty(F) # If the face set is empty, return all emtpy outputs
        Fc = F
        Vc = Vector{Point{ND,TV}}(undef,0)
        indFix = Vector{Int64}(undef,0)
    else # Faces not empty to check which indices are used and shorten V if needed        
        indUsed = elements2indices(F) # Indices used
        Vc = V[indUsed] # Remove unused points    
        indFix = zeros(Int64,length(V))
        indFix[indUsed] .= 1:length(indUsed)
        Fc = [(eltype(F))(indFix[f]) for f in F] # Fix indices in F         
    end
    return Fc, Vc, indFix
end

function trisurfslice(F::Vector{TriangleFace{TF}},V::Vector{Point{ND,TV}},n = Vec{3, Float64}(0.0,0.0,1.0), p = mean(V,dims=1); snapTolerance = 0.0, output_type=:full) where TF<:Integer where ND where TV<:Real 

    if !in(output_type,(:full,:above,:below))
        throw(ArgumentError("Invalid output_type :$output_type provided, use :full,:above, or :below"))
    end
    intersectFunc(v1,v2,d,n) = v1 .- d/dot(n,v2.-v1) .* (v2.-v1)
    
    # Compute dot product with normal of slicing plane
    d = map(v-> dot(n,v.-p),V)
    if snapTolerance != 0.0
        d[abs.(d).<snapTolerance] .= 0.0
    end
    LV = d.<0.0
    
    Fn =  Vector{TriangleFace{TF}}()
    Cn =  Vector{Int64}()
    Vn = deepcopy(V)
    D = Dict{Vector{Int64},Int64}() # For pointing from edge to intersection point index
    for f in F
        lf = LV[f]
        
        if any(lf) # Some or all below
            if all(lf) # All below
                if output_type == :full || output_type == :below            
                    push!(Fn,f)
                    push!(Cn,-2)
                end                 
            else # Some below -> cut
                nBelow = sum(lf) # Number of nodes below
                if isone(nBelow) # 2-above, 1 below 
                    indP = f[wrapindex(findfirst(lf) .+ (0:2),3)]
                    
                    e1 = sort(indP[[1,2]])
                    if !haskey(D,e1)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[2]],d[indP[1]],n))
                        D[e1] = length(Vn)
                    end
                    
                    e2 = sort(indP[[1,3]])
                    if !haskey(D,e2)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[3]],d[indP[1]],n))
                        D[e2] = length(Vn)
                    end

                    if output_type == :above || output_type == :full       
                        push!(Fn,TriangleFace{TF}(D[e1],indP[2],indP[3]))
                        push!(Fn,TriangleFace{TF}(D[e1],indP[3],D[e2]))
                        push!(Cn,1)
                        push!(Cn,1)
                    end
                    
                    if output_type == :below || output_type == :full
                        push!(Fn,TriangleFace{TF}(indP[1],D[e1],D[e2]))
                        push!(Cn,-1)                        
                    end

                else # 1-above, 2 below
                    indP = f[wrapindex(findfirst(.!lf) .+ (0:2),3)]

                    e1 = sort(indP[[1,2]])
                    if !haskey(D,e1)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[2]],d[indP[1]],n))
                        D[e1] = length(Vn)
                    end                    

                    e2 = sort(indP[[1,3]])
                    if !haskey(D,e2)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[3]],d[indP[1]],n))
                        D[e2] = length(Vn)
                    end
                    
                    if output_type == :below || output_type == :full                        
                        push!(Fn,TriangleFace{TF}(D[e1],indP[2],indP[3]))
                        push!(Fn,TriangleFace{TF}(D[e1],indP[3],D[e2]))
                        push!(Cn,-1)
                        push!(Cn,-1)
                    end

                    if output_type == :above || output_type == :full
                        push!(Fn,TriangleFace{TF}(indP[1],D[e1],D[e2]))
                        push!(Cn,1)                        
                    end
                end
            end
        else # Not any below -> all above
            if output_type == :full || output_type == :above            
                push!(Fn,f)
                push!(Cn,2)
            end    
        end
    end    
    Fn,Vn = remove_unused_vertices(Fn,Vn)
    return Fn,Vn,Cn
end

function count_edge_face(F,E_uni=nothing,indReverse=nothing)::Vector{Int64}
    if isnothing(E_uni) || isnothing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    con_F2E = con_face_edge(F,E_uni,indReverse)
    
    C = zeros(Int64,length(E_uni))
    for i_f in eachindex(F)
        for i in con_F2E[i_f]
            C[i]+=1
        end
    end 
    return C
end

function boundaryedges(F::Vector{NgonFace{N,TF}}) where N where TF <: Integer
    E = meshedges(F)
    Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
    count_E2F = count_edge_face(F,Eu,indReverse)
    return Eu[isone.(count_E2F)]
end

function edges2curve(Eb::Vector{LineFace{TF}}) where TF <: Integer
    # TO DO: 
    # Handle while loop safety/breaking
    # Cope with non-ordered meshes (normals not coherent) 

    con_E2E = con_edge_edge(Eb) # Edge-edge connectivity
    seen = fill(false,length(Eb)) # Bool to keep track of visited points
    i = 1 # Start with first edge
    ind = [Eb[i][1]] # Add first edge point and grow this list
    while !all(seen) # loop until all edges have been visited        
        push!(ind,Eb[i][2]) # Add edge end point (start is already in list)
        seen[i] = true # Lable current edge as visited       
        e_ind = con_E2E[i] # Indices for connected edges
        if Eb[e_ind[1]][1]==ind[end] #Check if 1st point of 1st edge equals end
            i = e_ind[1]
        elseif length(e_ind)>1 && Eb[e_ind[2]][1]==ind[end] #Check if 1st point of 2nd edge equals end
            i = e_ind[2]
        end
    end
    return ind
end

"""
    pointspacingmean(V::Vector{Point3{Float64}})
    pointspacingmean(F::Array{NgonFace{N, Int64}, 1},V::Vector{Point3{Float64}}) where N

The `pointspacingmean` function computes the mean spacing between points. The 
input can be just the coordinate set `V`, a vector of Point3 
points, or also a set of edges `E` or faces `F`. If only `V` is provided it is 
assumed that `V` represents an ordered set of "adjacent" points, e.g. as for a 
curve. If a vector of edges `E` or a vector of faces `F is also provided, then 
the average edge length is computed. If instead a set of faces `F` is provided 
then edges are first computed after which the mean edge spacing is return. 
"""
function pointspacingmean(V::Vector{Point{ND,TV}}) where ND where TV <: Real
    # Equivalent to:  mean(norm.(diff(V,dims=1)))
    p = 0.0
    n = length(V)
    @inbounds for i in 1:(n-1)
        p += norm(V[i]-V[i+1])/(n-1)
    end
    return p
end

function pointspacingmean(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real
    if isa(F,Vector{LineFace{Int64}})
        E = F
    else
        E = meshedges(F; unique_only=true)
    end
    p = 0.0
    n = length(E)
    @inbounds for i in 1:n
        p += norm(V[E[i][1]]-V[E[i][2]])/n
    end
    return p
end

function pointspacingmean(M::GeometryBasics.Mesh)    
    return pointspacingmean(faces(M),coordinates(M))
end


function extrudecurve(V1::Vector{Point{ND,TV}},d; s=1, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing,close_loop=false,face_type=:quad) where ND where TV<:Real
    # Derive num_steps from curve point spacing if missing    
    if isnothing(num_steps)
        num_steps = ceil(Int64,d/pointspacingmean(V1))
        if face_type==:tri
            num_steps = num_steps + Int64(iseven(num_steps)) # Force uneven
        end
    end
    # Create offset point depending on direction of extrude
    if isone(s) # Allong n from V1
        p = d.*n
    elseif isone(-s) # Against n from V1
        p = -d.*n
        V1 = reverse(V1)
    elseif iszero(s) # Extrude both ways from V1
        p = d.*n
        V1 = [(eltype(V1))(v.-p./2) for v in V1] #Shift V1 in negative direction
    end
    V2 = [(eltype(V1))(v.+p) for v in V1]  
    return loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type=face_type)
end


function meshgroup(F; con_type = :v)

    if con_type == :v # Group based on vertex connectivity 
        con_F2F = con_face_face_v(F)
    elseif con_type == :e # Group based on edge connectivity
        # EDGE-VERTEX connectivity
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

        # FACE-EDGE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse)    

        # EDGE-FACE connectivity
        con_E2F = con_edge_face(F,E_uni)

        # FACE-FACE connectivity
        con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)
    else 
        throw(ArgumentError("Wrong con_type :$con_type used, use :v or :e"))
    end

    if all(isempty.(con_F2F)) # Completely disconnected face set (e.g. raw STL import)
        C = collect(1:length(F))
    else
        C = fill(0,length(F))
        i = 1
        c = 1
        C[i] = c 
        seen = Set{Int64}(1)
        while length(seen)<length(F)
            np = length(seen)
            con_f2f = con_F2F[i]
            if !isempty(con_f2f)
                ind_F = reduce(vcat,con_f2f)
                i = Vector{Int64}()
                for ii in  ind_F
                    if !in(ii,seen)
                        push!(seen,ii)
                        C[ii] = c
                        push!(i,ii)
                    end
                end
            else
                if !in(i,seen)
                    push!(seen,i)
                    C[i] = c 
                end
            end
            if np == length(seen) # Group full
                c += 1 # Increment group counter
                i = findfirst(iszero.(C))
            end
        end
    end
    return C
end

function distmarch(F,V::Vector{Point{ND,TV}},indStart; d=nothing, dd=nothing, dist_tol=1e-3,con_V2V=nothing,l=nothing) where ND where TV<:Real

    # Get vertex-vertex connectivity
    if isnothing(con_V2V)
        con_V2V = con_vertex_vertex_f(F,V) 
    end

    # Compute "Laplacian umbrella" distances
    if isnothing(dd)
        dd = Dict{Vector{Int64},Float64}()  
        for i in eachindex(V)
            for ii in con_V2V[i]
                k = sort([i,ii])
                if !haskey(dd,k)
                    dd[sort(k)] = norm(V[i]-V[ii])
                end 
            end
        end
    end

    # Get/allocate distance vector
    if isnothing(d)
        d = fill(Inf,length(V))
    end

    if isnothing(l)
        l = fill(0,length(V))
    end

    # Set start distances to zero 
    is_isolated =  isempty.(con_V2V)
    d[is_isolated] .= NaN # Set isolated (non-connected) points to NaN
    d[indStart] .= 0.0
    l[indStart] .= 1:length(indStart)
    
    notGrowing = false
    dist_sum_previous = -1.0 # Set negative initially 
    count_inf_previous = length(d)-length(indStart) # number of Inf values currently
    while true                          
        for i in eachindex(V) # For each point            
            for ii in con_V2V[i] # Check umbrella neighbourhood
                # Get closest point and distance from umbrella
                minVal,minInd = findmin([d[ii],dd[sort([i,ii])]+d[i]])            
                if minInd==2
                    d[ii] = minVal # Distance                          
                    l[ii] = l[i] # Index
                end
            end            
        end
        bool_inf = isinf.(d) # Booling to check number of points left at Inf
        count_inf = count(bool_inf)
        if count_inf == count_inf_previous #!any(isinf.(d)) # Start checking once all are no longer Inf
            dist_sum = sum(d[.!is_isolated .&& .!bool_inf])
            if notGrowing # If we were here before
                if abs(dist_sum-dist_sum_previous)<dist_tol                                        
                    break                    
                end
            end
            notGrowing = true # Flip to denote we've been here           
            dist_sum_previous = dist_sum # Now start computing the sum to check convergence
        end
        count_inf_previous = count_inf
    end
    d[isinf.(d)] .= NaN # Change Inf to NaN
    return d,dd,l
end

# function distseedpoints(F,V,numPoints; ind=[1],dist_tol=1e-3)
    
#     con_V2V = con_vertex_vertex_f(F,V) 
#     d,dd,l = distmarch(F,V,ind; dist_tol=dist_tol,con_V2V=con_V2V)

#     if numPoints>1
#         @showprogress 1 "<distseedpoints>: Seeding points..." for q in 2:numPoints            
#             push!(ind,findmax(d)[2])
#             d,dd,l = distmarch(F,V,ind; dist_tol=dist_tol, dd=dd,d=d,con_V2V=con_V2V,l=l)        
#         end
#     end
#     return ind,d,l
# end


"""
    ray_triangle_intersect(F::Vector{TriangleFace{Int64}},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))
    ray_triangle_intersect(f::TriangleFace{Int64},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))

# Description 
This function can compute triangle-ray or triangle-line intersections through 
the use of the "Möller-Trumbore triangle-ray intersection algorithm" [1]. The 
required inputs are as follows: 

`F` an single face or a vector of faces, e.g. `Vector{TriangleFace{Int64}}`
`V` The triangle vertices as a vector of points, i.e. `Vector{Point{3, Float64}}`
`ray_vector` The ray vector which can be `Vector{Point{3, Float64}}` or `Vec3{Float64}`

The following optional input parameters can be provided: 
`rayType = :ray` (default) or `:line`. This defines wether the vector is treated as a ray (extends indefinately) or as a line (finite length)
`triSide = 1` (default) or `0` or `-1`. 
When `triSide=1` only the inward intersections are considered, e.g. when the ray or line enters the shape (ray/line is pointing against face normal)
When `triSide=-1` only the outward intersections are considered, e.g. when the ray or line exits the shape (ray/line is pointing allong face normal)
When `triSide=0` both inward and outward intersections are considered.
`tolEps = eps(Float64)` (default) 

# References 
1. [Möller, Tomas; Trumbore, Ben (1997). _Fast, Minimum Storage Ray-Triangle Intersection_. Journal of Graphics Tools. 2: 21-28. doi: 10.1080/10867651.1997.10487468.](https://doi.org/10.1080/10867651.1997.10487468)
"""
function ray_triangle_intersect(F::Vector{TriangleFace{TF}},V::Vector{Point{ND,TV}},ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64)) where TF <: Integer where ND where TV<:Real
    P = Vector{Point{ND,TV}}()
    indIntersect = Vector{Int64}()
    for qf in eachindex(F)
        p = ray_triangle_intersect(F[qf],V,ray_origin,ray_vector; rayType = rayType, triSide = triSide, tolEps = tolEps)        
        if !any(isnan.(p))
            push!(P,p)
            push!(indIntersect,qf)
        end
    end
    return P,indIntersect
end

function ray_triangle_intersect(f::TriangleFace{Int64},V::Vector{Point{ND,TV}},ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64)) where ND where TV<:Real

    # Edge vectors
    P1 = V[f[1]] # First corner point
    vec_edge_1 = V[f[2]].-P1 # Edge vector 1-2
    vec_edge_2 = V[f[3]].-P1 # Edge vector 1-3

    # Determine if ray/lines is capable of intersecting based on direction
    ray_cross_e2 = cross(ray_vector,vec_edge_2) 
    det_vec = dot(vec_edge_1,ray_cross_e2)  # Determinant det([-n' P21' P31'])
    if triSide == 1 # Pointing at face normals
        boolDet = det_vec>tolEps
    elseif triSide == 0 # Both ways
        boolDet = abs(det_vec)>tolEps
    elseif triSide == -1 # Pointing allong face normals
        boolDet = det_vec<tolEps
    end

    p = Point{ND,TV}(NaN,NaN,NaN)
    if boolDet        
        s = ray_origin.-P1
        u = dot(s,ray_cross_e2)/det_vec    
        if u >= 0 && u <= 1 # On triangle according to u            
            s_cross_e1 = cross(s,vec_edge_1)
            v = dot(ray_vector,s_cross_e1)/det_vec
            if v >= 0 && (u+v) <= 1 # On triangle according to both u and v
                # Allong ray/line coordinates i.e. intersection is at ray_origin + t.*ray_vector 
                t = dot(vec_edge_2,s_cross_e1)/det_vec                      
                if rayType == :ray || (rayType == :line && t>=0 && t<=1.0)                                                   
                    p = ray_origin .+ t.*ray_vector # same as: push!(P, P1 .+ u.*P21 .+ v.*P31)            
                end
            end
        end    
    end    
    return p 
end

"""
    mesh_curvature_polynomial(F::Vector{TriangleFace{Int64}},V::Vector{Point3{Float64}})
    mesh_curvature_polynomial(M::GeometryBasics.Mesh)

# Description
This function computes the mesh curvature at each vertex for the input mesh 
defined by the face `F` and the vertices `V`. A local polynomial is fitted to 
each point's "Laplacian umbrella" (point neighbourhood), and the curvature of 
this fitted form is derived. Instead of the mesh faces and vertices one may 
instead specify the `GeometryBasics.Mesh` `M` as the input. 

The reference below [1] provides more detail on the algorithm. In addition, this 
implementation was created with the help of [this helpful document](https://github.com/alecjacobson/geometry-processing-curvature/blob/master/README.md), 
which features a nice overview of the theory/steps involved in this algorithm. 

# References 
1. [F. Cazals and M. Pouget, _Estimating differential quantities using polynomial fitting of osculating jets_, Computer Aided Geometric Design, vol. 22, no. 2, pp. 121-146, Feb. 2005, doi: 10.1016/j.cagd.2004.09.004](https://doi.org/10.1016/j.cagd.2004.09.004)
"""
function mesh_curvature_polynomial(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real
    # Get the unique mesh edges
    E_uni = meshedges(F;unique_only=true) 

    # Get the vertex-to-vertex connectivity, i.e. the "Laplacian umbrellas"
    con_V2V = con_vertex_vertex(E_uni,V)

    
    NV = vertexnormal(F,V) # The vertex normal directions
    nz = Vec{3,Float64}(0.0,0.0,1.0) # A z-axis vector

    K1 = Vector{Float64}(undef,length(V)) # Allocate first principal curvature
    K2 = Vector{Float64}(undef,length(V)) # Allocate second principal curvature
    U1 = Vector{Vec3{Float64}}(undef,length(V)) # Allocate first principal curvature vector
    U2 = Vector{Vec3{Float64}}(undef,length(V)) # Allocate second principal curvature vector
    for q in eachindex(V)
        n = NV[q] # The current vertex normal
        Q = rotation_between(n,nz) # The rotation between the current normal and the z-axis
        ind = con_V2V[q] # The indices for the current Laplacian umbrella       
        vr = [Q*(v-V[q]) for v in V[ind]] # Rotate point set to a 2D problem

        # Set up polynomial fit
        T = Matrix{Float64}(undef,(length(ind),5))
        w = Matrix{Float64}(undef,(length(ind),1))
        for i = 1:length(ind)
            T[i,:] = [vr[i][1],vr[i][2],vr[i][1]^2,vr[i][1]*vr[i][2],vr[i][2]^2]
            w[i] = vr[i][3]
        end     
        a = T\w  # x = A\B solves the system of linear equations A*x = B

        E = 1.0 + a[1]^2
        F = a[1]*a[2]
        G = 1.0 + a[2]^2
        d = sqrt(a[1]^2+1.0+a[2]^2)
        e = (2.0*a[3]) / d
        f =       a[4] / d
        g = (2.0*a[5]) / d

        S = -[e f; f g] * inv([E F; F G])
        k,u = eigen(S) # Eigen decomposition to get first/second eigenvalue and vectors
        
        # Store derived quantities
        K1[q] = k[2]
        K2[q] = k[1]
        U1[q] = Q'*Vec3{Float64}(u[1,2],u[2,2],0.0)
        U2[q] = Q'*Vec3{Float64}(u[1,1],u[2,1],0.0)    
    end

    H = 0.5 * (K1.+K2) # Mean curvature
    G = K1.*K2 # Gaussian curvature

    return K1,K2,U1,U2,H,G
end
    
function mesh_curvature_polynomial(M::GeometryBasics.Mesh) 
        return mesh_curvature_polynomial(faces(M),coordinates(M))
end

"""
    separate_vertices(F::Array{NgonFace{N, Int64}, 1},V::Array{Point{M, T}, 1}) where N where M where T<:Real
    separate_vertices(M::GeometryBasics.Mesh)

This function takes the input mesh defined by the faces `F` and vertices `V` and
separates any shared vertices. It does this by giving each face its own set of 
unshared vertices. Note that any unused points are not returned in the output 
point array `Vn`. 
"""
function separate_vertices(F::Vector{NgonFace{N, TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real
    Vn = Vector{eltype(V)}()
    Fn = deepcopy(F)
    c = 0 
    for q in eachindex(F)
        f = F[q]
        m = length(f)
        Fn[q] = c .+ (1:m)
        c += m
        append!(Vn,V[f])
    end
    return Fn,Vn
end

function separate_vertices(M::GeometryBasics.Mesh)
    F = faces(M)
    V = coordinates(M)
    Fn,Vn = separate_vertices(F,V)    
    return GeometryBasics.Mesh(Vn,Fn)
end

"""
    curve_length(V::Vector{Point{ND,TV}}; close_loop=false) where ND where TV<:Real

This function computes the stepwise length of the input curve defined by the ND 
points in `V`. The output is a vector containining the distance for each point, 
and the total length therefore the last entry. 

If the optional parameter `close_loop` is set to `true` then it is assumed that
the curve should be seen as closed, i.e. the last entry is for returning to the 
start point from the last point in `V`. 
"""
function curve_length(V::Vector{Point{ND,TV}}; close_loop=false) where ND where TV<:Real
    if close_loop 
        return pushfirst!(cumsum(push!(norm.(diff(V)),norm(V[1]-V[end]))),0.0) # Along curve distance from start-to-start
    else
        return pushfirst!(cumsum(norm.(diff(V))),0.0) # Along curve distance from start-to-end        
    end 
end


"""
    evenly_sample(V::Vector{Point{ND,TV}}, n::Int64; rtol = 1e-8, niter = 1) where ND where TV<:Real

Evenly samples curves. 

# Description

This function aims to evenly resample the input curve defined by the ND points 
`V` using `n` points. The function returns the resampled points as well as the 
spline interpolator `S` used. The output points can also be retriebed by using: 
`S.(range(0.0, 1.0, n))`. 
Note that the even sampling is defined in terms of the curve length for a 4th 
order natural B-spline that interpolates the input data. Hence if significant 
curvature exists for the B-spline between two adjacent data points then the 
spacing between points in the output may be non-uniform (despite the allong 
B-spline distance being uniform). 
"""
function evenly_sample(V::Vector{Point{ND,TV}}, n::Int64; rtol = 1e-8, niter = 1) where ND where TV<:Real
    m = length(V)
    LL = curve_length(V) # Initialise as along curve (multi-linear) distance
    LL ./= last(LL) # Normalise
    S = BSplineKit.interpolate(LL, V, BSplineOrder(4), BSplineKit.Natural()) # Create interpolator
    for _ in 1:niter
        dS = Derivative() * S  # spline derivative
        L = similar(LL) # Initialise spline length vector
        L[1] = 0
        for i in 2:m
            # Compute length of segment [i-1, i]
            segment_length, _ = quadgk(LL[i-1], LL[i]; rtol) do t
                norm(dS(t))  # integrate |S'(t)| in segment [i, i + 1]
            end        
            L[i] = L[i - 1] + segment_length
        end
        L ./= last(L) # Normalise to 0-1 range
        S = BSplineKit.interpolate(L, V, BSplineOrder(4), BSplineKit.Natural()) # Create interpolator
    end
    l = range(0.0, 1.0, n) #Even allong curve distance 
    return S.(l), S # Evaluate interpolator at even distance increments
end

function evenly_space(V::Vector{Point{ND,TV}}, pointSpacing::T; rtol = 1e-8, niter = 1) where ND where TV<:Real where T<:Real    
    n = ceil(Int64,maximum(curve_length(V))/pointSpacing)
    V,_ = evenly_sample(V,n; rtol=rtol, niter=niter)
    return V
end

"""
    invert_faces(F::Vector{NgonFace{N, TF}, 1}) where N where TF<:Integer

Flips face orientations.

# Description

This function inverts the faces in `F`, such that the face normal will be 
flipped, by reversing the node order for each face. 
"""
function invert_faces(F::Vector{NgonFace{N, TF}}) where N where TF<:Integer     
    return map(f-> reverse(f),F) # [NgonFace{N, Int64}(reverse(f)) for f in F]    
end


"""
    R = kabsch_rot(V1::Array{Point{N, T}, 1},V2::Array{Point{N, TT}, 1}) where N where T<:Real where TT<:Real

# Description  
Computes the rotation tensor `R` to rotate the points in `V1` to best match the 
points in `V2`. 

# Reference
[Wolfgang Kabsch, _A solution for the best rotation to relate two sets of vectors_, Acta Crystallographica Section A, vol. 32, no. 5, pp. 922-923, 1976, doi: 10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873) 
[https://en.wikipedia.org/wiki/Kabsch_algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm) 
"""
function kabsch_rot(V1::Vector{Point{ND,TV1}},V2::Vector{Point{ND,TV2}}) where ND where TV1<:Real where TV2<:Real
    # Centre on means 
    V1 = V1.-mean(V1)
    V2 = V2.-mean(V2)
        
    # Compute A matrix
    A = zeros(Float64,3,3)
    for i = 1:3    
        for j= 1:3
            for q in eachindex(V1)
                @inbounds A[i,j] += V1[q][i]*V2[q][j] 
            end
        end
    end

    # Kabsch algorithm
    U, _, V = svd(A)    
    d = det(U)*det(V) #sign(det(U'*V))
    D = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    D[3,3] = d 
    return Rotations.RotMatrix3{Float64}(V*D*U')
end


"""
    F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist = 0, close_loop=true)   

# Description
This function implements swept lofting. The start curve `V1` is pulled allong the 
guide curve `Vc` while also gradually (linearly) morphing into the end curve 
`V2`. 
The optional parameter `face_type` (default :quad) defines the type of mesh 
faces uses. The same face types as `loftlinear` and `extrudecurve` are supported, 
i.e. `:quad`, `:tri_slash`, `tri`, or `quad2tri`. 
The optional parameter `num_twist` (default is 0) can be used to add an integer 
number (negative or positive) of full twists to the loft. 
Finally the optional parameter `close_loop` (default is `true`) determines if the
section curves are deemed closed or open ended. 
"""
function sweeploft(Vc::Vector{Point{ND,TV}},V1::Vector{Point{ND,TV}},V2::Vector{Point{ND,TV}}; face_type=:quad, num_twist = 0, close_loop=true) where ND where TV<:Real   
    np = length(V1) # Number of section points
    nc = length(Vc) # Number of curve steps

    # Centre start/end sections arond means (should be curve start/end points)
    V1b = [v.-Vc[1] for v in V1] 
    V2b = [v.-Vc[end] for v in V2] 

    # Determine rotation between sections 
    n3_1 = facenormal([collect(1:length(V1))],V1)[1] # normalizevector(Vc[2]-Vc[1])
    n1_1 = normalizevector(V1b[1]) 
    n2_1 = normalizevector(cross(n1_1,n3_1))
    n1_1 = normalizevector(cross(n3_1,n2_1))
    S1p = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    n3_2 = facenormal([collect(1:length(V2))],V2)[1] # normalizevector(Vc[end]-Vc[end-1])
    n1_2 = normalizevector(V2b[1]) 
    n2_2 = normalizevector(cross(n1_2,n3_2))
    n1_2 = normalizevector(cross(n3_2,n2_2))        
    S2p = mapreduce(permutedims,vcat,[n1_2,n2_2,n3_2])

    Q12 = RotMatrix3{Float64}(S1p\S2p) # nearest_rotation(S1p\S2p)
    
    # Rotate V2b to start orientation
    V2b = [Q12*v for v in V2b] 

    # Linearly loft "alligned" sections to create draft intermediate sections
    F,V = loftlinear(V1b,V2b;num_steps=nc,close_loop=close_loop,face_type=face_type)
    
    # Rotating and positioning all sections 
    n1p = n1_1
    for q = 1:nc # For all curve points 
        
        ind = (1:np) .+ (q-1)*np # Indices of the points to map
        if q==1 # Just take the start section when we are at the start
            V[ind] = V1
        elseif q == nc # Just take the end section when we are at the start
            V[ind] = V2            
        else # Rotate and position intermediate section 
            n3 = normalizevector(normalizevector(Vc[q]-Vc[q-1]) .+ normalizevector(Vc[q+1]-Vc[q]))                        
            n2 = normalizevector(cross(n1p,n3))
            n1 = normalizevector(cross(n3,n2))
            S2 = mapreduce(permutedims,vcat,[n1,n2,n3])

            Q12 = RotMatrix3{Float64}(S2\S1p)
            V[ind] = [(Q12*v).+Vc[q] for v in V[ind]]
            n1p = n1
        end   

        if q == nc-1 # Once here, second to last, a potential rotational mismatch needs to be resolve
            Q_fix = RotMatrix3{Float64}(S2\S2p) # Rotation between last and second to last
            t_a = Rotations.params(AngleAxis(Q_fix)) # Angle/axis representation
            if dot(t_a[2:end],n3)<0                 
                β_fix = t_a[1]
            else # Flip angle sign in this case 
                β_fix = -t_a[1]
            end

            if β_fix>pi # Use shorter negative direction instead
                β_fix = -(2.0*pi-β_fix)
            end
        
            # Distribute the fix as around curve rotations for each intermediate step            
            β_range = range(0,β_fix+(2.0*π*num_twist),nc) 
            for q = 2:nc-1
                ind = (1:np) .+ (q-1)*np
                ns = normalizevector(normalizevector(Vc[q]-Vc[q-1]) .+ normalizevector(Vc[q+1]-Vc[q]))
                Q = AngleAxis(β_range[q], ns[1], ns[2], ns[3])                  
                V[ind] = [(Q*(v-Vc[q])).+Vc[q] for v in V[ind]]
            end
        end
    end      
    return F,V
end

"""
    revolvecurve(Vc::Vector{Point{ND,TV}},θ=2.0*pi; s=0, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing,close_loop=false,face_type=:quad)  where ND where TV<:Real   

Revolves curves to build surfaces 

# Description

This function rotates the curve `Vc` by the angle `θ`, in the direction `s`,
around the vector `n`, to build the output mesh defined by the faces `F` and 
vertices `V`. 
"""
function revolvecurve(Vc::Vector{Point{ND,TV}},θ=2.0*pi; s=0, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing,close_loop=false,face_type=:quad)  where ND where TV<:Real   
    # Compute num_steps from curve point spacing
    if isnothing(num_steps)
        rMax = 0.0
        for v in Vc
            rNow = dot(normalizevector(cross(cross(n,v),n)),v)
            if !isnan(rNow)
                rMax = max(rMax,rNow)
            end
        end
        num_steps = ceil(Int64,(rMax*θ)/pointspacingmean(Vc))        
    end
    
    # Set up angle range
    if isone(s) # Positive direction
        θ_range = range(0,θ,num_steps) 
        Vc = reverse(Vc) # To keep normal direction consistent    
    elseif iszero(s) # Both positive and negative directions 
        θ_range = range(0,θ,num_steps) # Positive range 
        Q = AngleAxis(-θ/2, n[1], n[2], n[3])
        Vc = [Q*v for v in Vc] #Rotate in negative direction by half the angle
        Vc = reverse(Vc) # To keep normal direction consistent        
    elseif isone(-s) # Negative direction
        θ_range = range(0,-θ,num_steps)            
    end

    V = Vector{eltype(Vc)}()
    for θ in θ_range
        Q = AngleAxis(θ, n[1], n[2], n[3])
        Vn = [Q*v for v in Vc]
        append!(V,Vn)
    end    
   
    return loftpoints2surf(V,num_steps;close_loop=close_loop,face_type=face_type)
end


"""
    batman(n::Int64)
# Description
The `batman` function creates points on the curve for the Batman logo. The curve
is useful for testing surface meshing algorithms since it contains sharp 
transitions and pointy features. The user requests `n` points on the curve. The
default forces exactly `n` points which may result in an assymetric curve. To 
instead force symmetry the user can set the optional parameter `symmetric=true`. 
In this case the output will be symmetric allong the y-axis, however the number
of points on the curve may have increased (if the input `n` is not even). The
second optional input is the direction of the curve, i.e. if it is clockwise, 
`dir=:cw` or anti-clockwise `dir=:acw`. 
The implementation is based on a "parameterised Batman equation" [1](https://www.desmos.com/calculator/ajnzwedvql).
The following modifications where made, the curve is here centered around 
[0,0,0], scaled to be 2 in width, resampled evenly, and the default curve 
direction is anti-clockwise. 

# References 
1. https://www.desmos.com/calculator/ajnzwedvql
"""
function batman(n::Int64; symmetric = false, dir=:acw)
    tt = range(8,0.0,n)

    x = [( 0.3*t + 0.2*abs(t-1.0) + 2.2*abs(t-2.0) - 2.7*abs(t-3.0) -3.0*abs(t-5.0) + 3*abs(t-7.0) 
        + 5.0*sin(π/4.0*(abs(t-3.0)-abs(t-4.0)+1.0))
        + 5.0/4.0 * (abs(t-4.0)-abs(t-5.0)-1.0)^3.0 
        - 5.3*cos((π/2.0 + asin(47/53))*(abs(t-7.0) 
        - abs(t-8.0) - 1.0)/2.0) + 2.8 ) for t in tt] 

    y = [( 3.0/2.0*abs(t-1.0) - 3.0/2.0*abs(t-2.0) - 29.0/4.0*abs(t-4.0) + 29.0/4.0*abs(t-5.0) 
        + 7.0/16.0*(abs(t-2.0) - abs(t-3.0) -1.0)^4.0 
        + 4.5*sin( π/4.0*(abs(t-3.0) - abs(t-4.0) - 1.0) )
        - ((3.0/5.0*sqrt(2)) * abs(abs(t-5.0)-abs(t-7.0))^(5.0/2.0))
        + 6.4 * sin( ( (π/2.0 + asin(47/53)) * ((abs(t-7.0) 
        - abs(t-8.0) + 1.0)/2.0) ) + asin(56/64) ) + 4.95 ) for t in tt]   

    if symmetric == true    
        if iseven(n)
            m = round(Int64,n/2+1)
        else
            m = ceil(Int64,n/2)+1 
        end
        V,_ = evenly_sample([Point{3, Float64}(x[i]/22,y[i]/22,0.0) for i in eachindex(x)],m)
        V2 = [Point{3, Float64}(-V[i][1],V[i][2],0.0) for i in length(V)-1:-1:2]
        append!(V,V2)    
    else
        V = [Point{3, Float64}(x[i]/22,y[i]/22,0.0) for i in eachindex(x)]    
        V2 = [Point{3, Float64}(-V[i][1],V[i][2],0.0) for i in length(V)-1:-1:2]
        append!(V,V2) 
        V,_ = evenly_sample(V,n)
    end
    if dir==:cw
        reverse!(V)   
        circshift!(V,1)     
    end
    V = V.-mean(V,dims=1)
    return V
end

"""
    tridisc(r=1.0,n=0)

# Description

Generates the faces `F` and vertices `V` for a triangulated disc (circle). The 
algorithm starts with a triangulated hexagon (obtained if `n=0`) and uses 
iterative subtriangulation circular coordinate correction to obtain the final 
mesh. 
"""
function tridisc(r=1.0,n=0)
    # Create a triangulated hexagon
    V = circlepoints(r,6)
    push!(V,Point{3,Float64}(0.0,0.0,0.0))
    F = [TriangleFace{Int64}(1,2,7), TriangleFace{Int64}(2,3,7), TriangleFace{Int64}(3,4,7), TriangleFace{Int64}(4,5,7), TriangleFace{Int64}(5,6,7), TriangleFace{Int64}(6,1,7)]

    # Refine mesh n times
    if n>0
        @inbounds for _ = 1:n            
            F,V = subtri(F,V,1; method = :linear)
            indB = elements2indices(boundaryedges(F))
            V[indB] = r.*V[indB]./norm.(V[indB]) # Push to conform to circle 
        end
        indB = elements2indices(boundaryedges(F))
    end
    return F,V
end

# """
#     triplate(xSpan,ySpan,pointSpacing::T) where T <: Real    

# # Description
# Generates a triangulated mesh for a plate.
# """
# function triplate(xSpan,ySpan,pointSpacing::T) where T <: Real    
#     return gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = true, rectangular=true)
# end

"""
    regiontrimesh(VT,R,P)

# Description

Generates a multi-region triangle mesh for the input regions. The boundary 
curves for all regions are containedin the tuple `VT`. Each region to be meshed
is next defined using a tuple `R` containing indices into the curve typle `VT`. 
If an entry in `R` contains only one index then the entire curve domain is 
meshed. If `R` contains multiple indices then the first index is assumed to be 
for the outer boundary curve, while all subsequent indices are for boundaries 
defining holes in this region. 
"""
function regiontrimesh(VT,R,P)
    V = Vector{Point{3,Float64}}() # eltype(VT)()
    F = Vector{TriangleFace{Int64}}()
    C = Vector{Float64}()
    for q in eachindex(R)        
        r = R[q] # The curve indices for the current region       
        pointSpacing = P[q] # Region point spacing
        Vn = reduce(vcat,VT[r]) # The current curve point set 
        numCurvePoints = length(Vn) # Number of points in the current curve set

        # Creating constraints
        constrained_segments = Vector{Vector{Vector{Int64}}}()
        n = 1        
        for i in eachindex(r)
            m = length(VT[r[i]])            
            ind = append!(collect(n:(n-1)+m),n)
            if i>1 #&& i==length(r)
                reverse!(ind)
            end
            append!(constrained_segments,[[ind]])
            n += m
        end
        
        # Adding interior points 
        xSpan =[minimum([v[1] for v in Vn]),maximum([v[1] for v in Vn])]
        ySpan =[minimum([v[2] for v in Vn]),maximum([v[2] for v in Vn])]
        Vg = gridpoints_equilateral(xSpan,ySpan,pointSpacing)
        Vn = append!(Vn,Vg)
        
        # Initial triangulation 
        constrained_segments_ori = deepcopy(constrained_segments) # Clone since triangulate can add new constraint points
        TRn = triangulate(Vn; boundary_nodes=constrained_segments)
        Fn =  [TriangleFace{Int64}(tr) for tr in TRn.triangles]

        # Check if new boundary points were introduced and remove if needed 
        Eb = boundaryedges(Fn)
        indB = unique(reduce(vcat,Eb))
        indRemove = indB[indB.>numCurvePoints]        
        if !iszero(length(indRemove))            
            logicKeep = fill(true,length(Vn))
            logicKeep[indRemove] .= false
            indKeep = findall(logicKeep)
            indFix = zeros(Int64,length(Vn))
            indFix[indKeep] = 1:length(indKeep)
            Vn = Vn[indKeep]
            constrained_segments = [[indFix[c[1]]] for c in constrained_segments_ori] # Fix indices after point removal                
            TRn = triangulate(Vn; boundary_nodes=constrained_segments)
            Fn = [TriangleFace{Int64}(tr) for tr in TRn.triangles]
            Vn = TRn.points
        end

        # Remove unused points (e.g. outside region)
        Fn,Vn,indFix1 = remove_unused_vertices(Fn,Vn)
        constrained_segments = [[indFix1[c[1]]] for c in constrained_segments]  # Fix indices after point removal 

        # Remove 3 and 4 connected points
        E_uni = meshedges(Fn; unique_only=false)
        con_V2V = con_vertex_vertex(E_uni,Vn)
        nCon = map(length,con_V2V)
        Eb = boundaryedges(Fn)
        indB = unique(reduce(vcat,Eb))

        indLowCon = findall(nCon.<5)
        indRemove = indLowCon[ [!in(i,indB) for i in indLowCon] ]

        logicKeep = fill(true,length(Vn))
        logicKeep[indRemove] .= false
        indKeep = findall(logicKeep)
        indFix = zeros(Int64,length(Vn))
        indFix[indKeep] = 1:length(indKeep)

        Vn = Vn[indKeep] # Remove points 
        constrained_segments = [[indFix[c[1]]] for c in constrained_segments] # Fix indices after point removal 

        # Redo triangulation after points have been removed
        TRn = triangulate(Vn; boundary_nodes=constrained_segments)
        Fn = [TriangleFace{Int64}(tr) for tr in TRn.triangles]
        Vn = TRn.points
        Fn,Vn,indFix = remove_unused_vertices(Fn,Vn)

        # Smoothen mesh using Laplacian smoothing
        Eb = boundaryedges(Fn)
        indB = unique(reduce(vcat,Eb))
        n = 25
        λ = 0.5
        Vn = smoothmesh_laplacian(Fn,Vn, n, λ; constrained_points=indB)

        # Append to output 
        if !iszero(length(V))
            Fn = [f.+length(V) for f in Fn] # Shift indices of new faces
        end
        append!(V,Vn) # Append vertices 
        append!(F,Fn) # Append faces 
        append!(C,fill(q,length(Fn))) # Append color data 
    end
    F,V,_ = mergevertices(F,V) # merge shared vertices
    return F,V,C 
end

"""
    scalesimplex(F,V,s)
Scales faces (or general simplices) wrt their centre. 

# Description

This function scales each simplex (e.g. a face) wrt their centre (mean of 
coordinates). This function is useful in generating lattice structures from 
elements as well as to create visualisations whereby "looking into" the mesh is
needed. 
"""
function scalesimplex(F,V,s)
    F,V = separate_vertices(F,V) # separate_vertices
    for f in F 
        vm = mean(V[f],dims=1)
        V[f] = s.*(V[f].-vm).+vm        
    end
    return F,V
end

"""
    subcurve(V,n)

Adds points between each curve point.  

# Description

This function adds `n` points between each current point on the curve `V`.
"""
function subcurve(V,n)
    m = length(V)
    Vn = Vector{Point{3,Float64}}(undef,m+(m-1)*n)    
    for q = 1:lastindex(V)-1                
        vn = range(V[q],V[q+1],n+2)
        i1 = 1+(q-1)*(n+1)
        Vn[i1:i1+n] = vn[1:end-1]
    end
    Vn[end] = V[end]
    return Vn
end