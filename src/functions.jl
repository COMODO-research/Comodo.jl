# Define types
abstract type AbstractElement{N,T} <: StaticVector{N,T} end

GeometryBasics.@fixed_vector Element = AbstractElement

const Tet4{T} = Element{4,T} where T<:Integer
const Tet10{T} = Element{10,T} where T<:Integer
const Hex8{T} = Element{8,T} where T<:Integer
const Hex20{T} = Element{20,T} where T<:Integer
const Penta6{T} = Element{6,T} where T<:Integer
const Rhombicdodeca14{T} = Element{14,T} where T<:Integer
const Truncatedocta24{T} = Element{24,T} where T<:Integer


"""
    ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V, con_V2V_f, con_F2F_v)

# Description 

A struct featuring the connectivity data for a mesh.   
"""
struct ConnectivitySet{N}
    edge_vertex::Vector{LineFace{Int}}
    edge_face::Vector{Vector{Int}}
    edge_edge::Vector{Vector{Int}}
    face_vertex::Vector{NgonFace{N,Int}}
    face_edge::Vector{Vector{Int}}
    face_face::Vector{Vector{Int}}
    vertex_edge::Vector{Vector{Int}}
    vertex_face::Vector{Vector{Int}}
    vertex_vertex::Vector{Vector{Int}}
    vertex_vertex_f::Vector{Vector{Int}}
    face_face_v::Vector{Vector{Int}}
 end


"""
    comododir()

# Description 

This function simply returns the string for the Comodo path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function comododir()
    pkgdir(@__MODULE__)
end


"""
    slidercontrol(hSlider,ax)    

Adds arrow key control to sliders

# Description 

This function adds arrow key control to Makie sliders. The inputs are the 
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
function slidercontrol(hSlider::Slider,ax::Union{GLMakie.Axis, GLMakie.Axis3, GLMakie.Figure, GLMakie.LScene})    
    sliderRange = hSlider.range[] # Get slider range
    rangeLength = length(sliderRange) # Number of possible steps 
    sliderIndex = hSlider.selected_index[] # Current slider index
    GLMakie.on(events(ax).keyboardbutton) do event
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
    framerate = ceil(Int,length(stepRange)/duration)
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
    if isempty(F)
        return Set{eltype(eltype(F))}()
    else
        f = first(F)
        S = Set{eltype(f)}()
        for f in F 
            for j in f 
                push!(S, j)
            end
        end
        return S 
    end
end

struct Grid3D{T, X, Y, Z} <: AbstractVector{Point{3, T}}
    x::X
    y::Y
    z::Z
    function Grid3D(x::X, y::Y, z::Z) where {X, Y, Z}
        T = promote_type(eltype(x), eltype(y), eltype(z))
        return new{T, X, Y, Z}(x, y, z)
    end
end
Base.length(g::Grid3D) = length(g.x) * length(g.y) * length(g.z)
Base.size(g::Grid3D) = (length(g),)
function Base.getindex(g::Grid3D{T}, i::Int) where {T}
    nx, ny, nz = length(g.x), length(g.y), length(g.z)
    idx = CartesianIndices((1:nx, 1:ny, 1:nz))
    ix, iy, iz = Tuple(idx[i])
    return Point{3, T}(g.x[ix], g.y[iy], g.z[iz])
end

"""
    gridpoints(x::Vector{T}, y=x, z=x) where T<:Real

Returns 3D grids of points

# Description

The `gridpoints` function returns a vector of 3D points which span a grid in 3D 
space. Points are defined as per the input ranges or range vectors. The output 
point vector contains elements of the type `Point`. 
"""
function gridpoints(x, y=x, z=x) 
    return Grid3D(x, y, z)
end  

"""
    gridpoints_equilateral(xSpan::Union{Vector{TT},Tuple{TT,TT}},ySpan::Union{Vector{TT},Tuple{TT,TT}},pointSpacing::T; return_faces::Val{B1} = Val(false), rectangular::Val{B2}=Val(false), force_equilateral::Val{B3}=Val(false)) where {T<:Real, TT<:Real, B1, B2, B3}

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
function gridpoints_equilateral(xSpan::Union{Vector{TT},Tuple{TT,TT}},ySpan::Union{Vector{TT},Tuple{TT,TT}},pointSpacing::T; return_faces::Val{B1} = Val(false), rectangular::Val{B2}=Val(false), force_equilateral::Val{B3}=Val(false)) where {T<:Real, TT<:Real, B1, B2, B3}
    minX = minimum(xSpan)
    maxX = maximum(xSpan)
    minY = minimum(ySpan)
    maxY = maximum(ySpan)

    # Set up x-range
    wx = maxX - minX
    nx = ceil(Int,wx./pointSpacing)+1 # Number of points in the x-direction
    xRange = range(minX,maxX,nx) # The xRange
    pointSpacingReal_X = xRange[2]-xRange[1]
    
    # Set up y-range
    pointSpacing_Y = pointSpacingReal_X.*0.5*sqrt(3) # Point spacing adjusted for equilateral triangles
    if B3 # Perfect equilateral grid needed, does not conform to span        
        yRange = minY:pointSpacing_Y:maxY
    else # Approximate grid, conforms to span
        wy = maxY - minY        
        ny = ceil(Int,wy./pointSpacing_Y)+1 # Number of points in the y-direction
        yRange = range(minY,maxY,ny) # The yRange
    end
    
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
            if B2
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
    if B1
        plateElem=(length(xRange)-1,length(yRange)-1)
        F = Vector{TriangleFace{Int}}(undef,prod(plateElem)*2)
        num_x = length(xRange)
        ij2ind(i,j) = i + ((j-1)*num_x) # function to convert subscript to linear indices    
        c = 1
        for i = 1:plateElem[1]
            for j = 1:plateElem[2]      
                if iseven(j)  
                    F[c] = TriangleFace{Int}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1)])
                    c += 1
                    F[c] = TriangleFace{Int}([ij2ind(i+1,j+1),ij2ind(i,j+1),ij2ind(i,j)])
                    c += 1
                else
                    F[c] = TriangleFace{Int}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i,j+1)])
                    c += 1
                    F[c] = TriangleFace{Int}([ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)])
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
    if x isa AbstractRange
        xx = collect(x)
    else
        xx = deepcopy(x)
    end

    if y isa AbstractRange
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
        throw(ArgumentError("Invalid pad_data method provided, valid options are :linear, :constant, and :none"))
    end

    # Change behaviour depending on extrapolation method
    if extrapolate_method==:linear
        # Simple data based linear extrapolation 
        L = [xii<x[1] || xii>x[end] for xii in xi] # Boolean for points to extrapolate for
        # TODO: Optimise to remove need for masks
        if any(L) # If any points outside of the range were encountered
            yi = similar(xi) # Initialise yi
            yi[L] .= lerp(xx,yy,xi[L]) # Linearly extrapolate outside of the input range
            yi[.!L] = interp_biharmonic(xx,yy,xi[.!L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic(xx,yy,xi)
        end
    elseif extrapolate_method==:constant
        # Simple constant extrapolation (last/first value is repeated indefinately)
        Ls = xi.<x[1] # Boolean for points preceeding the start
        Ll = xi.>x[end] # Boolean for points after the end
        # TODO: Optimise to remove need for masks
        if any(_x -> _x < x[1] || _x > x[end], xi) # If any points outside of the range were encountered
            yi = similar(xi) # Initialise yi 
            copyto!(yi, xi)
            replace!(_x -> _x < x[1] ? yy[1] : _x > x[end] ? yy[end] : _x, yi) 

            L = .!Ls .&& .!Ll # Boolean for points to interpolate using biharmonic interpolation
            yi[L] = interp_biharmonic(xx,yy,xi[L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic(xx,yy,xi)
        end
    elseif extrapolate_method==:biharmonic
        # Allow extrapolation as per the biharmonic function
        yi = interp_biharmonic(xx,yy,xi) 
    else
        error("InvalidParameter: Invalid extrapolate_method method provided, valid options are :linear, :constant, and :biharmonic")
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
    replace!(x -> isnan(x) ? zero(x) : x, g) # Replace NaN entries by zeros
    W = g \ y # Weights  

    D = dist(xi,x) # Distance between points in X and XI

    G = (D.^2).*(log.(D).-1.0) # Green's function.
    replace!(x -> isnan(x) ? zero(x) : x, G) # Replace NaN entries by zeros
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
        throw(ArgumentError("n is too low. Request at least two data points"))
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
data `y` at the specified site `xi`. 
"""
function lerp(x, y, xi::T) where {T<:Real}
    if length(x) != length(y)
        throw(DimensionMismatch("x and y must be the same length"))
    end
    j = findfirst(>(xi), x)
    if isnothing(j)
        j = length(x)
        i = j-1
    elseif isone(j)
        i = 1
        j = 2
    else
        i = j-1
    end
    w = abs(x[j]-x[i])
    t = (xi-x[i])/w
    yi = (1.0-t)*y[i] + t*y[j]
    return yi
end
lerp(x, y, xi::AbstractVector) = lerp.((x,), (y,), xi)
  
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
    return pairwise(Euclidean(), V1, V2)
end

function dist(V1::Vector{T},v2::T) where T <: AbstractVector
    return pairwise(Euclidean(), V1, (v2,))
end

function dist(v1::T,V2::Vector{T}) where T <: AbstractVector
    return pairwise(Euclidean(), (v1,), V2)
end

"""
    mindist(V1,V2; getIndex=Val(false), skipSelf = false )

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
function mindist(V1,V2; getIndex::Val{B1}=Val(false), skipSelf = false ) where {B1}
    T = promote_type(eltype(eltype(V1)), eltype(eltype(V2)))
    D, d = similar(V1, T), similar(V2, T)
    if B1
        I = Vector{Int}(undef,length(V1))
    end
    for i in eachindex(V1)
        for j in eachindex(V2)
            if skipSelf && i==j
                d[j] = Inf
            else
                d[j] = euclidean(V1[i],V2[j]) # norm(V1[i]-V2[j]) 
            end       
        end
        if B1
            D[i], I[i] = findmin(d)
        else
            D[i] = minimum(d)
        end
    end
    if B1
        return D, I
    else
        return D
    end
end

"""
    occursonce(X::Union{Tuple{Vararg{T, N}}, Array{T, N}}; sort_entries=false) where T <: Any where N  

Checks if entries occur once

# Description
This function returns a Boolean vector where the entries are true where the 
corresponding entry in the input `X` is unique (occurs once). If the optional 
parameter `sort_entries` (default is false) is true, then each entry will be 
sorted, in this case and entry [3,1,2] is viewed as the same as [1,3,2] and 
[1,2,3] and so on. 
"""
function occursonce(X; sort_entries=false)  
    d = Dict{eltype(X), Int}()
    B = falses(length(X))
    for (i,x) in enumerate(X)    
        sort_entries && (x = _sort(x)) 
        if !haskey(d, x)            
            d[x] = i # index in dict            
            B[i] = true             
        else
            B[i] = false # Set current to false
            B[d[x]] = false # And also one found earlier           
        end
    end
    return B
end

@inline _sort(x) = length(x) > 1 ? sort(x) : x

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
function gunique(X; compute_unique::Val{CompUnique}=Val(true), compute_index::Val{CompIdx}=Val(false), compute_inverse::Val{CompInv}=Val(false), compute_counts::Val{CompCounts}=Val(false), sort_entries=false) where {CompUnique,CompIdx,CompInv,CompCounts}
    if CompUnique && !(CompIdx || CompInv || CompCounts)
        return sort_entries ? unique(_sort, X) : unique(X)
    else
        CompIdx && (unique_indices = Int[]; sizehint!(unique_indices, length(X)))
        CompInv && (inverse_indices = similar(X, Int))
        CompUnique && (unique_values = empty(X))
        CompCounts && (counts = Int[]; sizehint!(counts, length(X)))
        seen = Dict{eltype(X),Int}()
        sizehint!(seen, length(X))

        nunique = 0
        for (i, x) in enumerate(X)
            sort_entries && (x = _sort(x))
            if !haskey(seen, x)
                nunique += 1
                seen[x] = nunique
                CompUnique && (push!(unique_values, x))
                CompIdx && (push!(unique_indices, i))
                CompInv && (inverse_indices[i] = nunique)
                CompCounts && (push!(counts, 1))
            elseif CompInv || CompCounts
                idx = seen[x]
                CompInv && (inverse_indices[i] = idx)
                CompCounts && (counts[idx] += 1)
            end
        end

        ret = ()
        CompUnique && (ret = (unique_values,))
        CompIdx && (ret = (ret..., unique_indices))
        CompInv && (ret = (ret..., inverse_indices))
        CompCounts && (ret = (ret..., counts))

        return ret
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
    virtualFaceIndices = sub2ind(n.*ones(Int,length(F[1])),sort.(F))    
    _, ind1, ind2 = unique_dict(virtualFaceIndices) 

    return F[ind1], ind1, ind2
end


"""
    ind2sub(siz,ind)

Converts linear indices to subscript indices. 

# Description

Converts the linear indices in `ind`, for a matrix/array with size `siz`, to the 
equivalent subscript indices.  
"""
function ind2sub(siz::Union{Tuple{Vararg{Int, N}}, Array{Int, N}},ind::Union{Int,Tuple{Vararg{Int, M}}, Array{Int, M}}) where N where M
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
        A = Vector{Int}[]
    end
    return A
end

# ind2sub helper function to parse just a single linear index and produce a single subscript index set 
function ind2sub_(ind::Int,numDim::Int,k::Union{Int,Array{Int, N},Tuple{Vararg{Int, N}}}) where N
    a = Vector{Int}(undef,numDim) # Initialise a
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
function sub2ind(siz::Union{Tuple{Vararg{Int, N}}, Array{Int, N}},A::Union{Vector{Vector{Int}}, Array{NgonFace{M, Int}, 1}}) where N where M
    numDim = length(siz)
    k = cumprod([siz[i] for i in eachindex(siz)],dims=1)        
    ind = Vector{Int}(undef,length(A))
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

function sub2ind(siz::Union{Tuple{Vararg{Int, N}}, Vector{Int}},A::Vector{Int})  where N  
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

function sub2ind_(a::Union{Tuple{Vararg{Int, N}}, Vector{Int},NgonFace{M, Int}},numDim::Int,k::Union{Int,Vector{Int},Tuple{Vararg{Int, N}}})  where N where M
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
    nf = length(F)
    E = Vector{LineFace{T}}(undef,N*nf)        
    for (i,f) in enumerate(F) # Loop over each node/point for the current simplex                           
        E[i:nf:end] = [LineFace{T}(f[j],f[mod1(j+1,N)]) for j in 1:N]        
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
    F = Vector{TriangleFace{Int}}(undef,20)
    F[1 ] = TriangleFace{Int}(9,4,1)
    F[2 ] = TriangleFace{Int}(1,5,9)
    F[3 ] = TriangleFace{Int}(1,8,5)
    F[4 ] = TriangleFace{Int}(10,8,1)
    F[5 ] = TriangleFace{Int}(4,10,1)
    F[6 ] = TriangleFace{Int}(5,2,12)
    F[7 ] = TriangleFace{Int}(12,2,3)
    F[8 ] = TriangleFace{Int}(12,3,6)
    F[9 ] = TriangleFace{Int}(12,6,9)
    F[10] = TriangleFace{Int}(12,9,5)
    F[11] = TriangleFace{Int}(10,7,11)
    F[12] = TriangleFace{Int}(8,10,11)
    F[13] = TriangleFace{Int}(2,8,11)
    F[14] = TriangleFace{Int}(3,2,11)
    F[15] = TriangleFace{Int}(7,3,11)
    F[16] = TriangleFace{Int}(2,5,8)
    F[17] = TriangleFace{Int}(10,4,7)
    F[18] = TriangleFace{Int}(7,6,3)
    F[19] = TriangleFace{Int}(6,7,4)
    F[20] = TriangleFace{Int}(6,4,9)
    
    return F, V
end


"""
    octahedron(r=1.0)

Creates an octahedron mesh. 
    
# Description

Creates a GeometryBasics.Mesh for an octahedron with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function octahedron(r=1.0)
    # Create vertices    
    s = r/sqrt(2.0)
    V=Vector{Point{3, Float64}}(undef,6)
    V[1 ] = Point{3, Float64}(   -s,    -s,  0.0)
    V[2 ] = Point{3, Float64}(    s,    -s,  0.0)
    V[3 ] = Point{3, Float64}(    s,     s,  0.0)
    V[4 ] = Point{3, Float64}(   -s,     s,  0.0)
    V[5 ] = Point{3, Float64}(  0.0,   0.0,   -r)
    V[6 ] = Point{3, Float64}(  0.0,   0.0,    r)
    
    # Create faces
    F = Vector{TriangleFace{Int}}(undef,8)
    F[1 ] = TriangleFace{Int}(5,2,1)
    F[2 ] = TriangleFace{Int}(5,3,2)
    F[3 ] = TriangleFace{Int}(5,4,3)
    F[4 ] = TriangleFace{Int}(5,1,4)
    F[5 ] = TriangleFace{Int}(6,1,2)
    F[6 ] = TriangleFace{Int}(6,2,3)
    F[7 ] = TriangleFace{Int}(6,3,4)
    F[8 ] = TriangleFace{Int}(6,4,1)
    
    return F, V
end


"""
    dodecahedron(r=1.0)

Creates a dodecahedron mesh. 
    
# Description

Creates a GeometryBasics.Mesh for an dodecahedron with radius `r`. The default 
radius, when not supplied, is `1.0`. 
"""
function dodecahedron(r=1.0; h=1.0/Base.MathConstants.golden)
    # Setting up coordinate parameters
    # The h-factor is for the general pyritohedron with 1/ϕ being default for a
    # regular dodecahedron.

    s = r/sqrt(3.0) # Scaling factor
    t = s * (1.0 + h)
    w = s * (1.0 - h^2)

    # Create vertices    
    V = Vector{Point{3, Float64}}(undef,20)
    V[ 1] = Point{3, Float64}([ -s,  -s,  -s])
    V[ 2] = Point{3, Float64}([  s,  -s,  -s])
    V[ 3] = Point{3, Float64}([  s,   s,  -s])
    V[ 4] = Point{3, Float64}([ -s,   s,  -s])
    V[ 5] = Point{3, Float64}([ -s,  -s,   s])
    V[ 6] = Point{3, Float64}([  s,  -s,   s])
    V[ 7] = Point{3, Float64}([  s,   s,   s])
    V[ 8] = Point{3, Float64}([ -s,   s,   s])
    V[ 9] = Point{3, Float64}([ -w, 0.0,  -t])
    V[10] = Point{3, Float64}([  w, 0.0,  -t])
    V[11] = Point{3, Float64}([ -w, 0.0,   t])
    V[12] = Point{3, Float64}([  w, 0.0,   t])
    V[13] = Point{3, Float64}([0.0,  -t,  -w])
    V[14] = Point{3, Float64}([0.0,  -t,   w])
    V[15] = Point{3, Float64}([  t,   w, 0.0])
    V[16] = Point{3, Float64}([  t,  -w, 0.0])
    V[17] = Point{3, Float64}([0.0,   t,  -w])
    V[18] = Point{3, Float64}([0.0,   t,   w])
    V[19] = Point{3, Float64}([ -t,  -w, 0.0])
    V[20] = Point{3, Float64}([ -t,   w, 0.0])

    # Create faces
    F = Vector{NgonFace{5,Int}}(undef,12)
    F[ 1] = NgonFace{5,Int}([18,  8, 11, 12,  7])
    F[ 2] = NgonFace{5,Int}([12, 11,  5, 14,  6])
    F[ 3] = NgonFace{5,Int}([11,  8, 20, 19,  5])
    F[ 4] = NgonFace{5,Int}([20,  8, 18, 17,  4])
    F[ 5] = NgonFace{5,Int}([ 9,  1, 19, 20,  4])
    F[ 6] = NgonFace{5,Int}([19,  1, 13, 14,  5])
    F[ 7] = NgonFace{5,Int}([13,  1,  9, 10,  2])
    F[ 8] = NgonFace{5,Int}([13,  2, 16,  6, 14])
    F[ 9] = NgonFace{5,Int}([ 2, 10,  3, 15, 16])
    F[10] = NgonFace{5,Int}([ 7, 12,  6, 16, 15])
    F[11] = NgonFace{5,Int}([ 7, 15,  3, 17, 18])
    F[12] = NgonFace{5,Int}([10,  9,  4, 17,  3])
    
    return F, V
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
    F = Vector{QuadFace{Int}}(undef,6)
    F[1 ] = QuadFace{Int}(1,2,3,4)
    F[2 ] = QuadFace{Int}(8,7,6,5)
    F[3 ] = QuadFace{Int}(5,6,2,1)
    F[4 ] = QuadFace{Int}(6,7,3,2)    
    F[5 ] = QuadFace{Int}(7,8,4,3)    
    F[6 ] = QuadFace{Int}(8,5,1,4)    

    return F, V
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
    F = Vector{TriangleFace{Int}}(undef,4)
    F[1 ] = TriangleFace{Int}(1,2,3)
    F[2 ] = TriangleFace{Int}(4,2,1)
    F[3 ] = TriangleFace{Int}(4,3,2)
    F[4 ] = TriangleFace{Int}(4,1,3)
    
    return F,V
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
        F,V = tetrahedron(r)
    elseif n==2
        F,V = cube(r)
    elseif n==3
        F,V = octahedron(r)
    elseif n==4 
        F,V = icosahedron(r)
    elseif n==5
        F,V = dodecahedron(r)
    end
    return F,V
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
special integer type `OffsetInteger{-1, TF}` is converted to `Int`.  
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
        F = [LineFace{Int}(f) for f in FM]
    elseif m == 3 # Triangles
        F = [TriangleFace{Int}(f) for f in FM]
    elseif m == 4 # Quads
        F = [QuadFace{Int}(f) for f in FM]
    else # Other mesh type
        F = [NgonFace{m,Int}(f) for f in FM]
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

function topoints(VM::Array{Vec{N, T}, 1}) where T <: Real where N   
    return [Point{N, T}(v) for v in VM]
end

function topoints(VM::Vector{Vector{T}}) where T <: Real    
        m = length(VM[1])
        return [Point{m, T}(v) for v in VM]
end

function topoints(VM::Vector{Point{ND,TV}}) where ND where TV <: Real        
    return VM
end

#TEMP FIX
function topoints(VM)
    [Point{3,Float64}(v) for v in coordinates(VM)]
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

Returns the edge cross product, useful for nomal direction and area computations. 

    # Description

This function computes the so-called edge-cross-product for a input mesh that is
either defined by the faces `F` and vertices `V`. 
"""
function edgecrossproduct(F::Union{Vector{NgonFace{NF,TF}},Vector{Vector{TF}}},V::Vector{Point{ND,TV}}) where NF where TF<:Integer where ND where TV<:Real
    if isa(F,Vector{Vector{TF}})
        N = length(F[1])
    else
        N = NF
    end

    C = Vector{GeometryBasics.Vec{ND, TV}}(undef,length(F)) # Allocate array cross-product vectors      
    @inbounds for (q,f) in enumerate(F) # Loop over all faces
        c = zeros(Vec{3,Float64})
        @inbounds for q in 1:N # Loop from first to end-1            
            c  += cross(V[f[q]],V[f[mod1(q+1,N)]]) # Add next edge contribution          
        end
        C[q] = c./2 # Length = face area, direction is along normal vector
    end
    return C
end

function edgecrossproduct(f::Union{NgonFace{NF,TF},Vector{TF}},V::Vector{Point{ND,TV}}) where NF where TF<:Integer where ND where TV<:Real    
    if isa(f,Vector{TF})
        N = length(f)
    else
        N = NF
    end

    c = zeros(Vec{3,TV})
    @inbounds for q in 1:N # Loop from first to end-1            
        c += cross(V[f[q]],V[f[mod1(q+1,N)]]) # Add next edge contribution          
    end
    return c./2 # Length = face area, direction is along normal vector
end

"""
    facenormal(F,V)

Returns the normal directions for each face.

# Description

This function computes the per face normal directions for the input mesh defined 
either by the faces `F` and vertices `V` or by the GeometryBasics mesh `M`. 
"""
function facenormal(F::Union{Vector{NgonFace{NF,TF}},Vector{Vector{TF}}},V::Vector{Point{ND,TV}}) where NF where TF<:Integer where ND where TV<:Real
    C = edgecrossproduct(F,V)
    return C./norm.(C)
end

function facenormal(f::Union{NgonFace{NF,TF},Vector{TF}},V::Vector{Point{ND,TV}}) where NF where TF<:Integer where ND where TV<:Real    
    c = edgecrossproduct(f,V)
    return c./norm(c)
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
function subtri(F::Vector{NgonFace{3,TF}},V::Vector{Point{ND,TV}},n::Int; method = :linear, constrain_boundary=false) where TF<:Integer where ND where TV <: Real
    
    if iszero(n)
        return F,V
    elseif isone(n)
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        # Check for boundary edges        
        count_E2F = count_edge_face(F,Eu,indReverse)
        B_boundary = isone.(count_E2F)
        if any(B_boundary)
            treatBoundary = true
            Eb = Eu[B_boundary]
            indB = unique(reduce(vcat,Eb))
            con_V2V = con_vertex_vertex(Eb,V)
        else
            treatBoundary = false
        end

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
            for i in eachindex(Eu) # For each edge index       
                if treatBoundary && B_boundary[i]              
                    Vm[i] = 1/2 .*(V[Eu[i][1]] .+ V[Eu[i][2]]) # Normal mid-edge point
                else
                    F_touch = F[con_E2F[i]] 
                    indVerticesTouch = Vector{TF}() 
                    for f in F_touch        
                        b = f.!=Eu[i][1] .&& f.!=Eu[i][2]      
                        if any(b)  
                            append!(indVerticesTouch,f[b])           
                        end
                    end                   
                    Vm[i] = 3/8 .*(V[Eu[i][1]] .+ V[Eu[i][2]])  .+ 1/8 .* (V[indVerticesTouch[1]] .+ V[indVerticesTouch[2]])
                end
            end
    
            # Modified vertices for original vertices
            Vv = deepcopy(V)
            for i in eachindex(V)            
                if treatBoundary && in(i,indB)
                    if !constrain_boundary
                        Vv[i] = 6/8*V[i] + 1/8*(V[con_V2V[i][1]]+V[con_V2V[i][2]]) 
                    end
                else
                    B_vert_face = [any(f.==i) for f in F]
                    F_touch = F[B_vert_face] 
                    indVerticesTouch = Vector{TF}()
                    for f in F_touch                
                        indTouch = f[f.!=i]        
                        for i in indTouch 
                            if i ∉ indVerticesTouch 
                                push!(indVerticesTouch,i)
                            end
                        end
                    end
                    N = length(indVerticesTouch)                
                    v_sum = sum(V[indVerticesTouch],dims=1)[1]                
                    β = 1/N * (5/8-(3/8 +1/4*cos((2*π)/N))^2)        
                    Vv[i] = (1-N*β) .* V[i] .+ β*v_sum                
                end
            end    
            # Create complete point set
            Vn = [Vv;Vm] # Updated orignals and new "mid-edge-ish" points
        else
            throw(ArgumentError("Incorrect method :$method. Use :linear or :Loop"))
        end

        return Fn,Vn    

    elseif n>1
        for _ = 1:n
            F,V = subtri(F,V,1; method=method, constrain_boundary=constrain_boundary)
        end
        return F,V
    else
        throw(ArgumentError("n should be larger than or equal to 0"))
    end
end

"""
subquad(F::Vector{NgonFace{4,TF}},V::Vector{Point{ND,TV}},n::Int; method=:linear) where TF<:Integer where ND where TV <: Real
subquad(F::Vector{NgonFace{4,TF}},V::Vector{Point{ND,TV}},n::Int; method=:Catmull_Clark) where TF<:Integer where ND where TV <: Real

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
function subquad(F::Vector{NgonFace{4,TF}},V::Vector{Point{ND,TV}},n::Int; method=:linear, constrain_boundary=false) where TF<:Integer where ND where TV <: Real
    if iszero(n)
        return F,V
    elseif isone(n)

        # Get edges
        E = meshedges(F) # Non-unique edges
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        # Check for boundary edges        
        count_E2F = count_edge_face(F,Eu,indReverse)
        B_boundary = isone.(count_E2F)
        if any(B_boundary)
            treatBoundary = true
            Eb = Eu[B_boundary]
            indB = unique(reduce(vcat,Eb))
            con_V2V = con_vertex_vertex(Eb,V)
        else
            treatBoundary = false
        end

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
            Vf = simplexcenter(F,V) # Mid face points
            Ve_mid = simplexcenter(Eu,V) # Mid edge points

            # Edge points 
            Ve = Vector{Point{ND,TV}}(undef,length(Eu)) # Initialize edge points
            for q in eachindex(Eu)      
                if treatBoundary && B_boundary[q]              
                    Ve[q] = Ve_mid[q] # Normal mid-edge point
                else                   
                    Ve[q] = (mean(Vf[con_E2F[q]],dims=1)[1] .+ Ve_mid[q])./2.0
                end
            end

            # Vertex points 
            Vv = deepcopy(V) # Initialize vertex points
            for q in eachindex(V) # Loop over all vertices
                if treatBoundary && in(q,indB)
                    if !constrain_boundary
                        Vv[q] = 6/8*V[q] + 1/8*(V[con_V2V[q][1]]+V[con_V2V[q][2]]) 
                    end
                else
                    indF = con_V2F[q]
                    indE = con_V2E[q]
                    N = length(indF) # Number of faces (or edges) touching this vertex                    
                    Vv[q] = (mean(Vf[indF],dims=1)[1] .+ 2.0.*mean(Ve_mid[indE],dims=1)[1] .+ (N-3.0).*V[q])./N
                end
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
            for ii = 0:3
                Fn[q+ii*length(F)] = QuadFace{TF}([F[q][ii+1], con_F2E[q][ii+1]+nv, q+nv+ne, con_F2E[q][1+mod(3+ii,4)]+nv])                
            end            
        end

        return Fn,Vn
    elseif n>1
        for _ =1:n
            F,V = subquad(F,V,1;method=method, constrain_boundary=constrain_boundary)
        end
        return F,V
    else
        throw(ArgumentError("n should be larger than or equal to 0"))
    end
end


function pushtoradius_(v,r=1.0)    
    return v.*(r/norm(v)) # Push to sphere       
end


"""
    geosphere(n::Int,r::T; method=:linear) where T <: Real

Returns a geodesic sphere triangulation

# Description

This function returns a geodesic sphere triangulation based on the number of
refinement iterations `n` and the radius `r`. Geodesic spheres (aka Buckminster-Fuller
 spheres) are triangulations of a sphere that have near uniform edge lenghts. 
The algorithm starts with a regular icosahedron. Next this icosahedron is refined 
`n` times, while nodes are pushed to a sphere surface with radius `r` at each
iteration. Two methods are available, i.e. `:linear` (default) and `:Loop` 
(see also `subtri`). The former features simply linear splitting while the latter
features the Loop method which may produce a smoother result.
"""
function geosphere(n::Int,r::T; method=:linear) where T <: Real
    F,V = platonicsolid(4,r)    
    for _ = 1:n                
        # Set push start
        if method == :linear 
            s = length(V) # push from here as linear leaves original in place 
        elseif method == :Loop 
            s = 1 # push all as Loop alters coordinates
        end

        # Sub-divide sphere
        F,V = subtri(F,V,1; method=method)

        # Push altered/new points to sphere
        @inbounds for i in s:length(V) 
            V[i] = pushtoradius_(V[i],r) 
        end
    end
    return F,V
end

"""
    geosphere(n::Int,r::T) where T <: Real

Returns a geodesic sphere triangulation

# Description

This function returns a geodesic hemispherephere triangulation based on the number of
refinement iterations `n` and the radius `r`. Geodesic spheres (aka Buckminster-Fuller
 spheres) are triangulations of a sphere that have near uniform edge lenghts. 
The algorithm starts with a regular icosahedron that is adjusted to generate a half icosahedron. 
Next this icosahedron is refined  `n` times, while nodes are pushed to a sphere surface with radius `r` at each
iteration. 
"""
function hemisphere(n::Int,r::T; face_type=:tri, closed=false) where T <: Real    
    if n<0
        throw(ArgumentError("n should be >= 0"))
    end

    # Creating the faces and vertices of a full sphere mesh for "n=0"
    if face_type == :tri        
        F,V = geosphere(1,r) # A once refined icosahedron has an "equator" line

        # Rotating sphere so a vertex "points north" and equator is in xy-plane            
        Q = RotXYZ(0.0,atan(Base.MathConstants.golden),0.0) # Rotation matrix        
        V = [Point{3, Float64}(Q*v) for v in V] # Rotate vertices        
    elseif face_type == :quad        
        F,V = quadsphere(1,r) # A once refined cube has an "equator" line        
    else
        throw(ArgumentError("face_type should be :tri or :quad"))
    end

    # Now cut off bottom hemisphere    
    searchTol = r./1000.0 # Tolerance for cropping hemisphere (safe, somewhat arbitrary if much smaller then mesh edge lengths)
    F = [F[i] for i in findall(map(f-> mean([V[j][3] for j in f])>=-searchTol,F))] # Remove faces below equator
    F,V,_ = remove_unused_vertices(F,V) # Cleanup/remove unused vertices after faces were removed
    C = fill(1,length(F))

    if closed
        if face_type == :tri           
            Fb,Vb = tridisc(r,1; ngon=5, method=:linear, orientation=:down)
            Q = RotXYZ(0.0,0.0,pi/2.0) # Rotation matrix
            Vb = [Point{3, Float64}(Q*v) for v in Vb] # Rotate vertices   
        elseif face_type == :quad  
            Fb,Vb = quaddisc(r,0; method=:linear, orientation=:down)                   
        end
        
        # Add base
        indTopFaces = 1:length(F) # Indices of top
        append!(F, [f .+ length(V) for f in Fb])
        append!(V,Vb)
        append!(C,fill(2,length(Fb)))

        # Merge nodes
        V,_,indMap = mergevertices(V; pointSpacing = ((pi/2)*r)/(1+2^n))
        indexmap!(F,indMap)
    end
    
    # Refine sphere if needed 
    if n>0 # If larger then refining is needed       
        @inbounds for _ = 1:n # Refine once n times
            nv = length(V) # Number of points prior to refining
            nf = length(F)

            # Change refine behaviour based on face_type            
            if face_type == :tri                        
                F,V = subtri(F,V,1)            
            elseif face_type == :quad        
                F,V = subquad(F,V,1)                
            end
            C = repeat(C,outer=4)                        

            if closed
                indTopFaces = [j .+ i*nf  for i = 0:3 for j in indTopFaces]
                indPush = Vector{Int}()
                for f in F[indTopFaces]
                    for i in f
                        if i>nv
                            push!(indPush,i)
                        end
                    end
                end                 
            else
                indPush = nv+1:length(V)
            end

            # Now push newly introduced points to the sphere
            @inbounds for i in indPush
                V[i] = pushtoradius_(V[i],r) # Overwrite points
            end
        end        
    end

    return F,V,C
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
                
    
    E = Vector{Hex8{Int}}(undef,numElements) # Allocate elements

    @inbounds for q in 1:numElements
        ijk_1 = ind2sub(boxEl,ind1[q])    
        ijk_2 = ijk_1 .+ ijk_shift[2]
        ijk_3 = ijk_1 .+ ijk_shift[3]
        ijk_4 = ijk_1 .+ ijk_shift[4]
        ijk_5 = ijk_1 .+ ijk_shift[5]
        ijk_6 = ijk_1 .+ ijk_shift[6]
        ijk_7 = ijk_1 .+ ijk_shift[7]
        ijk_8 = ijk_1 .+ ijk_shift[8]
        
        E[q] = Hex8{Int}( sub2ind(boxNod,ijk_1),sub2ind(boxNod,ijk_2),sub2ind(boxNod,ijk_3),sub2ind(boxNod,ijk_4),
                            sub2ind(boxNod,ijk_5),sub2ind(boxNod,ijk_6),sub2ind(boxNod,ijk_7),sub2ind(boxNod,ijk_8) )
    end

    # Create vertices aka nodal coordinates
    indNodes = collect(Int,1:numNodes)
    IJK_nodes = ind2sub(boxNod,indNodes)

    V = convert(Vector{Point{3, Float64}},IJK_nodes)
    for q in eachindex(V)
        V[q] = ((V[q] .- Point{3, Float64}(1.0,1.0,1.0)).*(boxDim./boxEl)) .- Point{3, Float64}(boxDim/2.0)
    end

    # Create face sets from elements
    F = element2faces(E) 
    CF_type = repeat(collect(1:6),numElements) # Allocate face color/label data

    F_uni,indUni,c_uni = gunique(F,return_index=true, return_inverse=false,return_counts=true,sort_entries=true)
    Lb = isone.(c_uni)
    Fb = F_uni[Lb]
    CF_type_uni = CF_type[indUni]
    CFb_type = CF_type_uni[Lb]

    return E,V,F,Fb,CFb_type
end

"""
    con_face_edge(F,E_uni=nothing,indReverse=nothing)

Returns the edges connected to each face.

# Description

This function computes the face-edge connectivity. The input faces `F` (and 
optionally also the unique edges `E_uni` and reverse indices `indReverse` to map
to the non-unique edges, see also `gunique`) are used to create a list of edges 
connected to each face. If `F` contains N faces then the output contains N such 
lists. For triangles the output contains 3 edges per faces, for quads 4 per face
and so on.   
"""
function con_face_edge(F,E_uni=nothing,indReverse=nothing)
    if isnothing(E_uni) | isnothing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    return [Vector{Int}(a) for a in eachrow(reshape(indReverse,length(F),length(F[1])))] 
end


"""
con_edge_face(F,E_uni=nothing,indReverse=nothing)

Returns the faces connected to each edge.

# Description

This function computes the edge-face connectivity. The input faces `F` (and 
optionally also the unique edges `E_uni` and reverse indices `indReverse` to map
to the non-unique edges, see also `gunique`) are used to create a list of faces 
connected to each edges. If `E_uni` contains N edges then the output contains 
N such lists. For non-boundary edges each edge should connect to 2 faces. 
Boundary edges connect to just 1 face.  
"""
function con_edge_face(F,E_uni=nothing,indReverse=nothing)
    if isnothing(E_uni) || isnothing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    con_F2E = con_face_edge(F,E_uni,indReverse)
    
    con_E2F = [Vector{Int}() for _ in 1:length(E_uni)]
    for i_f in eachindex(F)
        for i in con_F2E[i_f]
            push!(con_E2F[i],i_f)
        end
    end 
    return con_E2F
end


"""
    con_face_face(F,E_uni=nothing,indReverse=nothing,con_E2F=nothing,con_F2E=nothing)

Returns the edge-connected faces for each face.

# Description

This function computes the face-face connectivity for each face. The input faces
`F` are used to create a list of faces connected to each face by a shared edge.
For non-boundary triangles for instance the output contains 3 edges per faces 
(which may be less for boundary triangles), and similarly non-boundary quads 
would each have 4 edge-connected faces. Additional optional inputs include: the 
unique edges `E_uni`, the reverse indices `indReverse` to map to the non-unique 
edges (see also `gunique`), as well as the edge-face `con_E2F` and face-edge 
`con_F2E` connectivity. These are all needed for computing the face-face 
connectivity and supplying them if already computed therefore save time.  
"""
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
        
        con_F2F = [Vector{Int}() for _ in 1:length(F)]        
        for i_f in eachindex(F)
            for i in reduce(vcat,con_E2F[con_F2E[i_f]])    
                if i!=i_f     
                    push!(con_F2F[i_f],i)
                end 
            end
        end        
        return con_F2F
    else # Just one face, so return empty
        return [Vector{Int}() ]
    end
end

"""
    con_face_face_v(F,V=nothing,con_V2F=nothing)

Returns the vertex-connected faces for each face.

# Description

This function computes the face-face connectivity for each face. The input faces
`F` are used to create a list of faces connected to each face by a shared vertex.
Additional optional inputs include: the vertices `V`, and the vertex-face 
connectivity `con_V2F`. In terms of vertices only the number of vertices, i.e. 
`length(V)` is neede, if `V` is not provided it is assumed that `length(V)` 
corresponds to the largest index in `F`. The vertex-face connectivity if not
supplied, will be computed by this function, hence computational time may be saved 
if it was already computed. 
"""
function con_face_face_v(F,V=nothing,con_V2F=nothing)
    if length(F)>1 # More than one face so compute connectivity
        if isnothing(con_V2F) 
            con_V2F = con_vertex_face(F,V)  # VERTEX-FACE connectivity
        end
        con_F2F = [Vector{Int}() for _ in 1:length(F)]
        for i_f in eachindex(F)
            for i in unique(reduce(vcat,con_V2F[F[i_f]]))    
                if i!=i_f     
                    push!(con_F2F[i_f],i)
                end 
            end
        end
        return con_F2F
    else # Just one face, so return empty
        return [Vector{Int}() ]
    end
end

"""
    con_vertex_simplex(F,V=nothing)

Returns how vertices connect to simplices

# Description

This function computes the vertex-simplex connectivity for each vertex. The input
simplices `F` are used to create a list of simplices connected to each vertex.
Additional optional inputs include: the vertices `V`.
In terms of vertices only the number of vertices, i.e. `length(V)` is needed, 
if `V` is not provided it is assumed that `length(V)` corresponds to the largest
index in `F`.  
"""
function con_vertex_simplex(F,V=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,F))
    else
        n = length(V)
    end
    con_V2F = [Vector{Int}() for _ in 1:n]
    for i_f in eachindex(F)
        for i in F[i_f]
            push!(con_V2F[i],i_f)
        end
    end
    return con_V2F
end

"""
    con_vertex_face(F,V=nothing)

Returns how vertices connect to faces

# Description
This function is an alias of `con_vertex_simplex`, and computes the vertex-face 
connectivity for each vertex. The input faces `F` are used to create a list of 
faces connected to each vertex. Additional optional inputs include: the vertices `V`.
In terms of vertices only the number of vertices, i.e. `length(V)` is needed, 
if `V` is not provided it is assumed that `length(V)` corresponds to the largest
index in `F`.  
"""
function con_vertex_face(F,V=nothing)
    return con_vertex_simplex(F,V)
end


"""
    con_vertex_edge(F,V=nothing)

Returns how vertices connect to edges

# Description
This function is an alias of `con_vertex_simplex`, and computes the vertex-edge 
connectivity for each vertex. The input edges `E` are used to create a list of 
edges connected to each vertex. Additional optional inputs include: the vertices `V`.
In terms of vertices only the number of vertices, i.e. `length(V)` is needed, 
if `V` is not provided it is assumed that `length(V)` corresponds to the largest
index in `E`.  
"""
function con_vertex_edge(E,V=nothing)
    return con_vertex_simplex(E,V)
end

"""
    con_edge_edge(E_uni,con_V2E=nothing)

Returns the vertex-connected edges for each edge.

# Description

This function computes the edge-edge connectivity for each edge. The input edges
`F` are used to create a list of edges connected to each edge by a shared vertex.
Additional optional inputs include: `con_V2E` (the vertex-edge connectivity), which
is instead computed when not provided.  
"""
function con_edge_edge(E_uni,con_V2E=nothing)
    if isnothing(con_V2E)
        con_V2E = con_vertex_edge(E_uni) 
    end    
    con_E2E = [Vector{Int}() for _ in 1:length(E_uni)]
    for i_e in eachindex(E_uni)
        for i in reduce(vcat,con_V2E[E_uni[i_e]])    
            if i!=i_e     
                push!(con_E2E[i_e],i)
            end 
        end
    end
    return con_E2E
end

"""
    con_vertex_vertex_f(F,V=nothing,con_V2F=nothing)

Returns the face-connected vertices for each vertex.

# Description

This function computes the vertex-vertex connectivity for each vertex using the 
vertex connected faces. The input faces `F` are used to create a list of vertices
connected to each vertex by a shared face. Additional optional inputs include: 
the vertices `V` and `con_V2F` (the vertex-face connectivity). In terms of vertices
only the number of vertices, i.e. `length(V)` is needed, if `V` is not provided 
it is assumed that `length(V)` corresponds to the largest index in `F`. The 
vertex-face connectivity `con_V2F` is needed, hence is computed when not provided.  
"""
function con_vertex_vertex_f(F,V=nothing,con_V2F=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,F))
    else 
        n = length(V)
    end
    
    if isnothing(con_V2F)
        con_V2F = con_vertex_face(F,V)
    end

    con_V2V = [Vector{Int}() for _ in 1:n]
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

"""
    con_vertex_vertex(E,V=nothing,con_V2E=nothing)

Returns the edge-connected vertices for each vertex.

# Description

This function computes the vertex-vertex connectivity for each vertex using the 
vertex connected edges. The input edges `E` are used to create a list of vertices
connected to each vertex by a shared edge. Additional optional inputs include: 
the vertices `V` and `con_V2E` (the vertex-edge connectivity). In terms of vertices
only the number of vertices, i.e. `length(V)` is needed, if `V` is not provided 
it is assumed that `length(V)` corresponds to the largest index in `E`. The 
vertex-edge connectivity `con_V2E` is needed, hence is computed when not provided.  
"""
function con_vertex_vertex(E,V=nothing,con_V2E=nothing)
    if isnothing(V)
        n = maximum(reduce(vcat,E))
    else 
        n = length(V)
    end
    if isnothing(con_V2E)
        con_V2E = con_vertex_edge(E,V)
    end

    con_V2V = [Vector{Int}() for _ in 1:n]
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

"""
    meshconnectivity(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real

Returns all mesh connectivity data

# Description
This function returns the `ConnectivitySet`, i.e. all mesh connectivity data for 
the input mesh defined by the faces `F` and the vertices `V`. 
The `ConnectivitySet` contains the following connectivity descriptions: 
* face-edge
* edge-face
* face-face
* face-face (wrt vertices)
* vertex-face
* vertex-edge
* edge-edge
* vertex-vertex
* vertex-vertex (wrt faces)
"""
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
    con_V2V_f = con_vertex_vertex_f(F,V,con_V2F)

    # FACE-FACE connectivity wrt vertices
    con_F2F_v = con_face_face_v(F,con_V2F)

    return ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V, con_V2V_f, con_F2F_v) 
end 

"""
    mergevertices(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; roundVertices = true, numDigitsMerge=nothing) where N where TF<:Integer where ND where TV<:Real

Merges points that coincide

# Description
This function take the faces `F` and vertices `V` and merges points that are sufficiently 
similar. Once points are merged the indices in `F` are corrected for the new reduced
point set. 
"""
function mergevertices(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; roundVertices = true, pointSpacing=nothing, numDigitsMerge=nothing) where N where TF<:Integer where ND where TV<:Real
    if isnothing(numDigitsMerge)         
        if isnothing(pointSpacing)
            pointSpacing = pointspacingmean(F,V)
        end
    end
    V,indUnique,indMap = mergevertices(V; roundVertices = roundVertices, pointSpacing=pointSpacing, numDigitsMerge=numDigitsMerge)
    indexmap!(F,indMap)
    return F,V,indUnique,indMap
end

function mergevertices(V::Vector{Point{ND,TV}}; roundVertices = true, pointSpacing=nothing, numDigitsMerge=nothing) where ND where TV<:Real
    if roundVertices
        if isnothing(numDigitsMerge) 
            # Compute numDigitsMerge from point spacing
            if isnothing(pointSpacing)
                # pointSpacing = norm(maxp(V)-minp(V))/length(V)
                 throw(ArgumentError("Specify either numDigitsMerge or pointSpacing"))
            else
                numDigitsMerge = 6-mag(pointSpacing)
            end            
        end        
        # Create rounded coordinates to help obtain unique set
        # Note -0.0+0.0 = 0.0 so addition of zero points helps avoid 0.0 and -0.0 being seen as unique
        VR = [round.(v,digits = numDigitsMerge)+Point{ND,TV}(0.0,0.0,0.0) for v in V]

        # Get unique indices and reverse for rounded vertices
        _,indUnique,indMap = gunique(VR; return_index=true, return_inverse=true,sort_entries=false)
        V = V[indUnique] # The unique node set
    else
        V,indUnique,indMap = gunique(V; return_index=true, return_inverse=true,sort_entries=false)
    end
    return V,indUnique,indMap
end

""" 
    indexmap!(F,indMap::Union{Vector{T},AbstractRange{T}}) where T<:Integer

Updates indices using map

# Description
This function assumes the first inputs `F` is a vector (or vector of vectors) 
containing integers such as indices. Next each i-th entry in `F`, denoted `f` is
updated using `F[i] = indMap[f]`. An example of the use of this function is 
in combination with mergevertices, which removed duplicated vertices. Once 
duplicates are removed the indices in a face array `F` may need to be updated. 
For instance nodes 2 and 5 in a list of 6 vertices is the same then 
`indMap = [1,2,3,4,2,5]`. Such that if F was `F=[[1,2,3],[4,5,6]]` it is mapped
to be: `F=[[1,2,3],[4,2,5]]`. 
"""
function indexmap!(F,indMap::Union{Vector{T},AbstractRange{T}}) where T<:Integer
    @inbounds for (i,f) in enumerate(F)
        F[i] = indMap[f] 
    end        
end

function indexmap(F,indMap::Union{Vector{T},AbstractRange{T}}) where T<:Integer
    FF = deepcopy(F)    
    indexmap!(FF,indMap)
    return FF    
end


"""
    mag(n::T) where T <: Real

Returns order of magnitude

# Description
This function returns the "order of magnitude" of the input number n in terms of
powers of 10. The order of magnitude is here defined as: 
` floor(Int,log10(abs(n)))`
Hence the order of the absolute number is returned and non-integer numbers are 
supported. 
Here is an example: 
`mag.([1, pi, 9.9, 10, 99, 100, 999.9, 1000])`
returns:  
     `[0,  0,   0, 1,   1,   2,     2,    3]`
"""
function mag(n::T) where T <: Real
    return floor(Int,log10(abs(n)))
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
function smoothmesh_laplacian(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}, n=1, λ=0.5; con_V2V=nothing, tolDist=nothing, constrained_points=nothing) where N where TF<:Integer where ND where TV<:Real
    
    if λ>1.0 || λ<0.0
        throw(ArgumentError("λ should be in the range 0-1"))
    end
    
    if λ>0.0
        if n==0
            return V
        elseif n>0
            indSmooth = elements2indices(F) # Indices of points involved in smoothing

            if maximum(indSmooth)>length(V) || minimum(indSmooth)<1
                throw(ErrorException("Out of range indices detected"))
            end
            
            # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if nothing
            if isnothing(con_V2V)
                E_uni = meshedges(F;unique_only=true)
                con_V2V = con_vertex_vertex(E_uni)
            end        
            c = 0
            while c<n             
                Vs = deepcopy(V)
                @inbounds for i in indSmooth # eachindex(V)
                    Vs[i] = (1.0-λ).*Vs[i] .+ λ*mean(V[con_V2V[i]]) # Linear blend between original and pure Laplacian
                end
                if !isnothing(constrained_points)
                    Vs[constrained_points] = V[constrained_points] # Put back constrained points
                end
                if !isnothing(tolDist) # Include tolerance based termination
                    d = 0.0
                    @inbounds for i in indSmooth # eachindex(V)
                        d+=sqrt(sum((V[i].-Vs[i]).^2)) # Sum of distances
                    end
                    if d<tolDist # Sum of distance smaller than tolerance?
                        break
                    end            
                end
                c+=1 
                V = Vs                
            end            
        else #n<0
            throw(ArgumentError("n should be greater or equal to 0"))
        end
    end
    return V
end


"""
    smoothmesh_hc(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}, n=1, α=0.1, β=0.5; con_V2V=nothing, tolDist=nothing, constrained_points=nothing) where N where TF<:Integer where ND where TV<:Real

# Description 

This function implements HC (Humphrey's Classes) smoothing [1]. This method uses
Laplacian like smoothing but aims to compensate for shrinkage/swelling by also 
"pushing back" towards the original coordinates. 

# References 
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
        indSmooth = elements2indices(F) # Indices of points involved in smoothing
        if maximum(indSmooth)>length(V) || minimum(indSmooth)<1
            throw(ErrorException("Out of range indices detected"))
        end

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
            @inbounds for i in indSmooth # eachindex(V)
                P[i] = mean(Q[con_V2V[i]]) # Laplacian 
                # Compute different vector between P and a point between original 
                # point and Q (which is P before laplacian)
                B[i] = P[i] .- (α.*V[i] .+ (1.0-α).*Q[i])
            end
            
            @inbounds for i in indSmooth # eachindex(V)      
                # Push points back based on blending between pure difference vector
                # B and the Laplacian mean of these      
                P[i] = P[i] .- (β.*B[i] .+ (1.0-β).* mean(B[con_V2V[i]]))
            end   
            
            if !isnothing(constrained_points)
                P[constrained_points] = V[constrained_points] # Put back constrained points
            end

            if !isnothing(tolDist) # Include tolerance based termination
                d = 0.0
                @inbounds for i in indSmooth #eachindex(V)
                    d+=sqrt(sum((P[i].-Q[i]).^2)) # Sum of distances
                end
                if d<tolDist # Sum of distance smaller than tolerance?
                    break
                end            
            end
            c+=1 
        end
        return P
    end
end


"""
    quadplate(plateDim,plateElem; orientation=:up)

Returns a quad mesh for a plate

# Description
This function creates a quadrilateral mesh (faces `F` and vertices `V`) for a 
plate. The dimensions in the x-, and y-direction are specified in the input vector 
`plateDim`, and the number of elements to use in each direction in the input 
vector `plateElem`. 
"""
function quadplate(plateDim,plateElem; orientation=:up)

    if  !in(orientation,(:up,:down))        
        throw(ArgumentError("Orientation not supported. Use :up or :down"))
    end 
    
    num_x = plateElem[1]+1
    num_y = plateElem[2]+1
    V = gridpoints(range(-plateDim[1]/2,plateDim[1]/2,num_x),range(-plateDim[2]/2,plateDim[2]/2,num_y),0.0)
    
    F = Vector{QuadFace{Int}}(undef,prod(plateElem))
    ij2ind(i,j) = i + ((j-1)*num_x) # function to convert subscript to linear indices    
    c = 1
    @inbounds for i = 1:plateElem[1]
        @inbounds for j = 1:plateElem[2]        
            F[c] = QuadFace{Int}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)])
            c += 1
        end
    end
    
    if orientation == :up
        return F, V
    else#if orientation == :down
        return invert_faces(F), V    
    end    
end

"""
    quadsphere(n::Int,r::T) where T <: Real

Returns a quadrangulated sphere

# Description
This function creates a quadrilateral mesh (faces `F` and vertices `V`) for a sphere
with a radius defined by the input `r`. The input `n` defines the density of sphere
mesh. The quad mesh is constructed using `subquad` subdivision of a regular cube, 
whereby `n` sets the number of splitting iterations to use. Using `n=0` therefore
returns a cube.
"""
function quadsphere(n::Int,r::T) where T <: Real
    F,V = platonicsolid(2,r)    
    if n > 0
        for _ in 1:n
            numPointsNow = length(V)
            F,V = subquad(F,V,1)
            @inbounds for i in numPointsNow+1:length(V)
                V[i] = pushtoradius_(V[i],r) # Overwrite points
            end
        end
    end
    return F,V
end

function simplex2vertexdata(F::Union{Vector{NgonFace{NF,TF}},Vector{Element{NE,TE}}},DF,V::Union{Nothing,Vector{Point{ND,TV}}}=nothing; con_V2F=nothing, weighting=:none) where NF where TF<:Integer where NE where TE<:Integer where ND where TV<:Real    
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
        if !isempty(con_V2F[q])
            if weighting==:none || length(con_V2F[q])==1 
                DV[q] = mean(T,DF[con_V2F[q]])
            elseif weighting==:area            
                a = A[con_V2F[q]]                
                DV[q] = sum(T,DF[con_V2F[q]].*a)./sum(a)
            end
        else
            DV[q] = T(NaN)
        end
    end
    return DV
end

function vertex2simplexdata(F,DV)
    T = eltype(DV) # Element type of data in DV
    DF =  (typeof(DV))(undef,length(F)) # Allocate data for F
    @inbounds for (i,f) in enumerate(F)        
        DF[i] = mean(view(DV,f)) # The mean of the vertex data for each entry in F
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

function circlepoints(r::T,n::Int; dir=:acw) where T <: Real
    return [Point{3, Float64}(r*cos(t),r*sin(t),0) for t in circlerange(n;dir=dir)]
end

function circlepoints(r::Tuple{T,T},n::Int; dir=:acw) where T <: Real
    return [Point{3, Float64}(r[1]*cos(t),r[2]*sin(t),0) for t in circlerange(n;dir=dir)]
end

function circlepoints(f::FunctionType,n::Int; dir=:acw) where {FunctionType <: Function}
    return [Point{3, Float64}(f(t)*cos(t),f(t)*sin(t),0) for t in circlerange(n;dir=dir)]
end

"""
    loftlinear(V1,V2;num_steps=2,close_loop=true,face_type=:tri)

Lofts surface between curves

# Description 

The `loftlinear` function spans a surface from input curve `V1` to curve `V2`. 
The surface is formed by "lerping" curves from `V1` to `V2` in `num_steps` 
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
function loftlinear(V1::Vector{Point{ND,TV}}, V2::Vector{Point{ND,TV}}; num_steps=nothing, close_loop=true, face_type=:quad) where ND where TV<:Real
    # Derive num_steps from distance and mean curve point spacing if missing    
    if isnothing(num_steps)
        d = mean([norm(V1[i]-V2[i]) for i in eachindex(V1)])
        dp = 0.5* (pointspacingmean(V1)+pointspacingmean(V2))
        num_steps = ceil(Int,d/dp)                
        if num_steps < 2
            num_steps = 2
        end
    end

    if num_steps < 2
        throw(ArgumentError("num_steps should be >=2"))
    end

    # Linearly blending points from first to last
    V = Vector{eltype(V1)}()
    for q in range(0,num_steps,num_steps)
        λ = q/num_steps
        Vn = (1.0-λ).*V1 .+ λ.* V2  
        append!(V,Vn)
    end   

    F = grid2surf(V,num_steps; face_type=face_type, periodicity=(close_loop,false))

    return F,V
end 

"""
    grid2surf(V::Vector{Point{ND,TV}},num_steps; face_type=:quad, periodicity=(false,false), tri_dir=2) where ND where TV<:Real

Creates surface from grid

# Description 
This function creates the faces `F` to form a closed surface for the input grid
defined by the vector of points `V`. The grid is assumed to be structured i.e. 
that it is composed of `num_steps` point sets which each have point-to-point 
correspondence. Such a grid can represent a set of corresponding curves, or a 
set of columns of points for a Cartesian grid for instance. 
Optional inputs include the following: 
* `face_type`, which can be `:quad` (default), `:tri`, `:tri_even`, `:backslash`, 
`:forwardslash`, or `:quad2tri`. 
* `periodicity`, which is of the type `Tuple{Bool, Bool}`. If 
`periodicity[1]==true` it is assumed the grid is periodic in the first direction
, i.e. that the grid should be closed in the first direction. If 
`periodicity[2]==true` it is assumed the grid is periodic in the second 
direction, i.e. that the grid should be closed in the second direction.
"""
function grid2surf(V::Vector{Point{ND,TV}},num_steps; face_type=:quad, periodicity=(false,false), tri_dir=1) where ND where TV<:Real
    
    # Get number of points in each offset curve
    nc = length(V)/num_steps # Number of points in curve
    if !isinteger(nc) || nc<1
        throw(ArgumentError("The length(V)/num_steps should produce an integer >1 but is instead $nc"))
    end
    nc = Int(nc)

    ij2ind(i,j) = i + nc*(j-1) # function to convert subscript to linear indices

    if periodicity[1]
        iEnd = nc
    else
        iEnd = nc-1
    end
    if periodicity[2]
        jEnd = num_steps
    else
        jEnd = num_steps-1
    end
    
    if face_type == :tri || face_type == :tri_even        

        if face_type == :tri
            f = 2.0 # divide by 2 to get half-spacing shift
            stepSize = 2 # Shift every second            
        else    
            f = 4.0 # divide by 4 to get quarter-spacing shift            
            stepSize = 1 # Shift every second            
        end

        if tri_dir == 1             
            if nc>3
                spline_order = 4            
            else
                spline_order = 2
            end            
            
            indEnd = num_steps 
            @inbounds for qq in stepSize:stepSize:indEnd
                indNow = (1:nc) .+ (qq-1)*nc

                L = curve_length(V[indNow])
                dL = diff(L) # Between point distances           
                if periodicity[1]
                    p = norm(V[indNow[1]]-V[indNow[end]]) 
                    push!(dL,p)                
                    if iseven(qq)
                        L_i = L .+ dL./f
                    else
                        L_i = L .- circshift(dL./f,1)
                    end
                else
                    if iseven(qq)
                        L_i = L[2:end-1] .+ dL[2:end]./f
                    else
                        L_i = L[2:end-1] .- dL[1:end-1]./f                
                    end
                end            
                if periodicity[1]                
                    bc = BSplineKit.Periodic(last(L) + p) # Use periodic bc for closed curves
                else
                    if spline_order == 4 
                        bc = BSplineKit.Natural() # Natural
                    else
                        bc = nothing
                    end
                end                  
                S = BSplineKit.interpolate(L, V[indNow], BSplineOrder(spline_order), bc) # Create interpolator

                # Even range for curve distance 
                if periodicity[1]
                    V[indNow] = S.(L_i)
                else
                    V[indNow[2:end-1]] = S.(L_i)
                end
            end
        elseif tri_dir == 2           
            if num_steps>3
                spline_order = 4            
            else
                spline_order = 2
            end            

            indEnd = nc 
            @inbounds for qq in stepSize:stepSize:indEnd
                indNow = qq:nc:nc*num_steps                        
                
                L = curve_length(V[indNow])
                dL = diff(L) # Between point distances           
                if periodicity[2]
                    p = norm(V[indNow[1]]-V[indNow[end]]) 
                    push!(dL,p)                
                    if iseven(qq)
                        L_i = L .+ dL./f
                    else
                        L_i = L .- circshift(dL./f,1)
                    end
                else
                    if iseven(qq)
                        L_i = L[2:end-1] .+ dL[2:end]./f
                    else
                        L_i = L[2:end-1] .- dL[1:end-1]./f                
                    end
                end            
                if periodicity[2]                
                    bc = BSplineKit.Periodic(last(L) + p) # Use periodic bc for closed curves
                else
                    if spline_order == 4                         
                        bc = BSplineKit.Natural() # Natural
                    else
                        bc = nothing
                    end
                end                                  
                S = BSplineKit.interpolate(L, V[indNow], BSplineOrder(spline_order), bc) # Create interpolator

                # Even range for curve distance 
                if periodicity[2]
                    V[indNow] = S.(L_i)
                else
                    V[indNow[2:end-1]] = S.(L_i)
                end
            end
        end
    end

    # Build faces
    i_f = 1
    if face_type == :quad || face_type == :quad2tri      
        F = Vector{QuadFace{Int}}(undef,iEnd*jEnd)             
        @inbounds for i = 1:iEnd
            @inbounds for j = 1:jEnd    
                F[i_f] = QuadFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j),ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)))
                i_f += 1
            end
        end
        if face_type == :quad2tri 
            F = quad2tri(F,V; convert_method = :angle)
        end
    elseif face_type == :tri || face_type == :tri_even
        F = Vector{TriangleFace{Int}}(undef,iEnd*jEnd*2)     
        if tri_dir == 2       
            @inbounds for j = 1:jEnd
                @inbounds for i = 1:iEnd    
                    if iseven(i)
                        F[i_f  ] = TriangleFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j), ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)) ) # 1-2-3
                        F[i_f+1] = TriangleFace{Int}( ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)), ij2ind(i,j) ) # 3-4-1
                    else
                        F[i_f  ] = TriangleFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j), ij2ind(i,mod1(j+1,num_steps)) ) # 1-2-4
                        F[i_f+1] = TriangleFace{Int}( ij2ind(mod1(i+1,nc),j), ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)) ) # 2-3-4
                    end
                    i_f +=2
                end
            end
        elseif tri_dir == 1
            @inbounds for i = 1:iEnd
                @inbounds for j = 1:jEnd    
                    if iseven(j)
                        F[i_f  ] = TriangleFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j), ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)) ) # 1-2-3
                        F[i_f+1] = TriangleFace{Int}( ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)), ij2ind(i,j) ) # 3-4-1
                    else
                        F[i_f  ] = TriangleFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j), ij2ind(i,mod1(j+1,num_steps)) ) # 1-2-4
                        F[i_f+1] = TriangleFace{Int}( ij2ind(mod1(i+1,nc),j), ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)) ) # 2-3-4
                    end
                    i_f +=2
                end
            end
        end

    elseif face_type == :backslash
        F = Vector{TriangleFace{Int}}(undef,iEnd*jEnd*2)     
        @inbounds for i = 1:iEnd
            @inbounds for j = 1:jEnd    
                F[i_f  ] = TriangleFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j), ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)) ) # 1-2-3
                F[i_f+1] = TriangleFace{Int}( ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)), ij2ind(i,j) ) # 3-4-1
                i_f +=2
            end
        end
    elseif face_type == :forwardslash
        F = Vector{TriangleFace{Int}}(undef,iEnd*jEnd*2)     
        @inbounds for i = 1:iEnd
            @inbounds for j = 1:jEnd    
                F[i_f  ] = TriangleFace{Int}( ij2ind(i,j), ij2ind(mod1(i+1,nc),j), ij2ind(i,mod1(j+1,num_steps)) ) # 1-2-4
                F[i_f+1] = TriangleFace{Int}( ij2ind(mod1(i+1,nc),j), ij2ind(mod1(i+1,nc),mod1(j+1,num_steps)), ij2ind(i,mod1(j+1,num_steps)) ) # 2-3-4
                i_f +=2
            end
        end    
    else
        throw(ArgumentError("Invalid face_type specified :$face_type, use :quad, :tri, :tri_even, :backslash, :forwardslash, or :quad2tri"))
    end

    return F
end


function dirplot(ax,V::Vector{Point{ND,TV1}},U::Union{Vector{Point{ND,TV2}},Vector{Vec{ND,TV2}}}; color=:black,linewidth=3,scaleval=1.0,style=:from) where ND where TV1 <: Real where TV2 <: Real
    E = [LineFace{Int}(i,i+length(V)) for i in 1:length(V)]

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

function dirplot(ax,V::Union{Point{ND,TV1},Vec{ND,TV1}},U::Union{Point{ND,TV2},Vec{ND,TV2}}; color=:black,linewidth=3,scaleval=1.0,style=:from) where ND where TV1 <: Real where TV2 <: Real
    
    E = [LineFace{Int}(1,2)]

    if style==:from        
        P = [V]
        push!(P,V.+(scaleval.*U))
    elseif style==:to
        P = [V.-(scaleval.*U)]
        push!(P,V)        
    elseif style==:through
        UU = (scaleval.*U)/2
        P = [V.-UU]
        push!(P,V.+UU)        
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

function edgeangles(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; deg=false) where N where TF<:Integer where ND where TV<:Real        
    A = Vector{GeometryBasics.Vec{N, Float64}}(undef,length(F))
    for (j,f) in enumerate(F)
        a = Vector{Float64}(undef,N)
        @inbounds for i in 1:N            
            ip1 = mod1(i+1,N)            
            ip2 = mod1(i-1,N)
            n1 = normalizevector(V[f[ip1]]-V[f[i]])
            n2 = normalizevector(V[f[ip2]]-V[f[i]])
            if deg 
                a[i] = acosd(clamp(dot(n1,n2),-1.0,1.0))
            else
                a[i] = acos(clamp(dot(n1,n2),-1.0,1.0))
            end
        end        
        A[j] = a        
    end
    return A
end

function quad2tri(F::Vector{QuadFace{TF}},V::Vector{Point{ND,TV}}; convert_method = :angle, eps_level=1e-9) where TF<:Integer where ND where TV<:Real
   
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
            if δab<(δaf-eps_level) # If significantly better (using eps_level avoids "noisy" numerical behaviour)
                ft = fb
            else
                ft = ff
            end    
        else
            throw(ArgumentError("Incorrect convert_method set $convert_method, use :forward, :backward, or :angle"))
        end
        push!(Ft,ft[1])
        push!(Ft,ft[2])
    end
    return Ft
end

function remove_unused_vertices(F,V::Vector{Point{ND,TV}}) where ND where TV<:Real
    if isempty(F) # If the face set is empty, return all emtpy outputs
        Fc = F
        Vc = Vector{Point{ND,TV}}()
        indFix = Vector{Int}()
    else # Faces not empty, so check which indices are used and shorten V if needed        
        indUsed = elements2indices(F) # Indices used
        Vc = V[indUsed] # Remove unused points    
        indFix = zeros(Int,length(V))
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
    Cn =  Vector{Int}()
    Vn = deepcopy(V)
    D = Dict{Vector{Int},Int}() # For pointing from edge to intersection point index
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
                    indP = f[mod1.(findfirst(lf) .+ (0:2),3)]
                    
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
                    indP = f[mod1.(findfirst(.!lf) .+ (0:2),3)]

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
    Fn,Vn,_ = remove_unused_vertices(Fn,Vn)
    return Fn,Vn,Cn
end

function count_edge_face(F,E_uni=nothing,indReverse=nothing)::Vector{Int}
    if isnothing(E_uni) || isnothing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    con_F2E = con_face_edge(F,E_uni,indReverse)
    
    C = zeros(Int,length(E_uni))
    for i_f in eachindex(F)
        for i in con_F2E[i_f]
            C[i]+=1
        end
    end 
    return C
end


"""
    boundaryedges(F::Vector{NgonFace{N,TF}}) where N where TF <: Integer

Returns boundary edges

# Description
This function returns the boundary edges `Eb` for the input mesh defined either 
by a vector of faces or edges. Both should be of the type NgonFace (e.g. 
LineFace for edges or TriangleFace, QuadFace, etc for faces). Boundary edges 
are those edges that occur only once in the mesh. 
"""
function boundaryedges(F::Vector{NgonFace{N,T}}) where N where T <: Integer
    if N == 2 # Input represents edges
        return F[occursonce(F; sort_entries=true)]
    else # Input assumed to represent faces
        E = meshedges(F)
        return E[occursonce(E; sort_entries=true)]
    end
end


"""
    boundaryfaces(F::Vector{NgonFace{N,TF}}) where N where TF <: Integer
    boundaryfaces(E::Vector{Element{N,T}}) where N where T <: Integer

Returns boundary faces

# Description
This function returns the boundary faces `Fb` for the input mesh defined either 
by a vector of elements or faces. 

Boundary faces are those faces that occur only once in the mesh. 
"""
function boundaryfaces(F::Vector{NgonFace{N,T}}) where N where T <: Integer
    return F[occursonce(F; sort_entries=true)]
end

function boundaryfaces(E::Vector{Element{N,T}}) where N where T <: Integer
    F = element2faces(E)
    return F[occursonce(F; sort_entries=true)]
end

"""
    boundaryfaceindices(F::Vector{NgonFace{N,T}}) where N where T <: Integer

Returns boundary face indices

# Description
This function returns the boundary face indices for the input mesh defined by 
a vector of faces. 

Boundary faces are those faces that occur only once in the mesh. 
"""
function boundaryfaceindices(F::Vector{NgonFace{N,T}}) where N where T <: Integer
    return findall(occursonce(F; sort_entries=true))
end

""" 
    edges2curve(Eb::Vector{LineFace{T}}) where T <: Integer

Converts boundary edges to a curve

# Description
This function takes a set of boundary edges `Eb`, which may not be ordered e.g.
consequtively, and returns an ordered set of indices `ind` defining a curve of 
consequtive points. The function returns empty output if the input edges do not 
for a proper curve. 
"""
function edges2curve(Eb::Vector{LineFace{T}}; remove_last = false) where T <: Integer
    # TO DO: 
    # Handle while loop safety/breaking
    # Cope with non-ordered meshes (normals not coherent) 

    if isempty(Eb)
        return Vector{T}[] # Return empty
    else
        numEdges = length(Eb)        
        con_V2E = con_vertex_edge(Eb) # Vertex-edge connectivity
        con_E2E = con_edge_edge(Eb,con_V2E) # Edge-edge connectivity
        numConnectedEdges = length.(con_V2E)
        indStartEnd = findall(numConnectedEdges .== 1)
        if length(indStartEnd) == 2 # Non-closed curve           
            if Eb[con_V2E[indStartEnd[1]][1]][1] == indStartEnd[1]
                i = con_V2E[indStartEnd[1]][1]
            else
                i = con_V2E[indStartEnd[2]][1]
            end
        elseif length(indStartEnd) == 0 # Consistent with closed loop
            i = 1
        else
            throw(ErrorException("Invalid edges. Edges may contain branches or multiple disconnected sets"))
        end

        seen = fill(false,numEdges) # Bool to keep track of visited points
        ind = [Eb[i][1]] # Add first edge point and grow this list
        while !all(seen) # loop until all edges have been visited        
            push!(ind,Eb[i][2]) # Add edge end point (start is already in list)
            seen[i] = true # Lable current edge as visited       
            e_ind = con_E2E[i] # Indices for connected edges
            if length(e_ind)>2 # Branch point detected
                throw(ErrorException("Invalid edges or branch point detected. Current edge is connected to more than two edges."))    
            else
                if Eb[e_ind[1]][1]==ind[end] #Check if 1st point of 1st edge equals end
                    i = e_ind[1]
                elseif length(e_ind)>1 && Eb[e_ind[2]][1]==ind[end] #Check if 1st point of 2nd edge equals end
                    i = e_ind[2]
                end
            end

            if seen[i] & !all(seen)
                throw(ErrorException("Invalid edges. Edges may contain multiple disconnected sets, such as multiple closed loops."))
            end            
        end
        if remove_last 
            return ind[1:end-1]
        else
            return ind
        end

    end
end


"""
    pointspacingmean(V::Vector{Point3{Float64}})
    pointspacingmean(F::Array{NgonFace{N, Int}, 1},V::Vector{Point3{Float64}}) where N

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
    if isa(F,Vector{LineFace{Int}})
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


"""
    extrudecurve(V1::Vector{Point{ND,TV}}; extent=1.0, direction=:positive, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing,close_loop=false,face_type=:quad) where ND where TV<:Real

Extrudes curves into surfaces
 
# Description
Extrudes (e.g. extends) the input curve defined by the points `V1` to a surface. 
The following input parameters are defined: 
* `extent<:Real` (default = 1.0) the length of the extrusion   
* `direction` is a symbol that is either `:positive` (default), `:negative`, or `:both`. 
* `n<:Vec{3, Float64}` The extrusion direction vector. The default is: Vec{3, Float64}(0.0,0.0,1.0). 
* `num_steps` (default is `nothing`) is the number of nodes in the extrude direction, the 
number of faces in the extrude direction is therefore `num_steps-1`. If not 
provided or `nothing` then the number of steps is derived so the point spacing in 
the extrusion direction matches the mean input curve point spacing. 
* `close_loop<:Bool` if true the curve is assumed to be closed, e.g. a curcle 
yields a closed cylinder
* `face_type` is a symbol that is either `:quad` (default), `tri`, `tri_slash`, 
or `quad2tri`. 
"""
function extrudecurve(V1::Vector{Point{ND,TV}}; extent=1.0, direction=:positive, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing,close_loop=false,face_type=:quad) where ND where TV<:Real

    if close_loop
        n1 = facenormal([collect(1:length(V1))],V1)[1] # Curve normal vector
        if dot(n,n1)<0 # Check if anti-aligned with extrusion direction
            V1 = circshift(reverse(V1),1) # Reverse curve order
        end
    end

    # Derive num_steps from curve point spacing if missing    
    if isnothing(num_steps)
        num_steps = ceil(Int,extent/pointspacingmean(V1))
        if face_type==:tri
            num_steps = num_steps + Int(iseven(num_steps)) # Force uneven
        end
        if num_steps < 2
            num_steps = 2
        end
    end

    # Check if num_steps is okay
    if num_steps<1
        throw(ArgumentError("num_steps=$num_steps is not valid. num_steps should be larger than 0."))
    end
    
    # Create offset point depending on direction of extrude
    if direction == :positive # Allong n from V1
        p = extent.*n
    elseif direction == :negative # Against n from V1
        p = -extent.*n
        if close_loop 
            circshift!(reverse!(V1),1)
        else
            reverse!(V1)
        end
    elseif direction == :both # Extrude both ways from V1
        p = extent.*n
        V1 = [(eltype(V1))(v.-p./2) for v in V1] #Shift V1 in negative direction
    else
        throw(ArgumentError("$direction is not a valid direction, Use :positive, :in, or :both.")) 
    end
    V2 = [(eltype(V1))(v.+p) for v in V1]  
    return loftlinear(V1,V2;num_steps=num_steps,close_loop=close_loop,face_type=face_type)
end

"""
   meshgroup(F; con_type = :v)

Groups connected mesh features
 
# Description
This function uses the connectivity `con_type` to create a group label `C` for
each entitiy in `F`. E.g. `C.==1` for the first group and `C.==2`` for the 
second and so on.     
"""
function meshgroup(F; con_type = :v, indStart=1, stop_at = nothing)

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
        C = fill(0,length(F)) # Initialise group label vector
        i = indStart
        c = 1
        C[i] = c 
        seen = Set{Int}(1)
        while length(seen)<length(F)
            np = length(seen)
            con_f2f = con_F2F[i]
            if !isempty(con_f2f)
                ind_F = reduce(vcat,con_f2f)
                i = Vector{Int}()
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
                if !isnothing(stop_at)                 
                    if c==stop_at
                        break
                    end                     
                end
                c += 1 # Increment group counter                
                i = findfirst(iszero.(C))                
                C[i] = c
            end
        end
    end
    return C
end

"""
    distmarch(F,V::Vector{Point{ND,TV}},indStart; d=nothing, dd=nothing, dist_tol=1e-3,con_V2V=nothing,l=nothing) where ND where TV<:Real

Compute on surface distance

# Description
This function computes allong mesh-edge distances for the points with the 
indices contained in `indStart`. 
"""
function distmarch(F,V::Vector{Point{ND,TV}},indStart; d=nothing, dd=nothing, dist_tol=1e-3,con_V2V=nothing,l=nothing) where ND where TV<:Real

    # Get vertex-vertex connectivity
    if isnothing(con_V2V)
        con_V2V = con_vertex_vertex_f(F,V) 
    end

    # Compute "Laplacian umbrella" distances
    if isnothing(dd)
        dd = Dict{Vector{Int},Float64}()  
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
    ray_triangle_intersect(F::Vector{TriangleFace{Int}},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))
    ray_triangle_intersect(f::TriangleFace{Int},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))

# Description 
This function can compute triangle-ray or triangle-line intersections through 
the use of the "Möller-Trumbore triangle-ray intersection algorithm" [1]. The 
required inputs are as follows: 

`F` an single face or a vector of faces, e.g. `Vector{TriangleFace{Int}}`
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
function ray_triangle_intersect(F::Vector{TriangleFace{TF}},V::Vector{Point{ND,TV1}},ray_origin::Union{Point{ND,TV2},Vec{ND,TV2}},ray_vector::Union{Point{ND,TV3},Vec{ND,TV3}}; rayType = :ray, triSide = 1, tolEps = eps(Float64)) where TF <: Integer where ND where TV1<:Real where TV2<:Real where TV3<:Real
    P = Vector{Point{ND,TV1}}()
    indIntersect = Vector{Int}()
    for (i,f) in enumerate(F)
        p = ray_triangle_intersect(f,V,ray_origin,ray_vector; rayType = rayType, triSide = triSide, tolEps = tolEps)        
        if !any(isnan.(p))
            push!(P,p)
            push!(indIntersect,i)
        end
    end
    return P,indIntersect
end

function ray_triangle_intersect(f::TriangleFace{Int},V::Vector{Point{ND,TV1}},ray_origin::Union{Point{ND,TV2},Vec{ND,TV2}},ray_vector::Union{Point{ND,TV3},Vec{ND,TV3}}; rayType = :ray, triSide = 1, tolEps = eps(Float64)) where ND where TV1<:Real where TV2<:Real where TV3<:Real

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

    p = Point{ND,TV1}(NaN,NaN,NaN)
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
    mesh_curvature_polynomial(F::Vector{TriangleFace{Int}},V::Vector{Point3{Float64}})
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
function mesh_curvature_polynomial(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}; growsteps = 2) where N where TF<:Integer where ND where TV<:Real

    if N>3 # e.g. Quads and up, which need face connectivity
        # Get the vertex-to-vertex connectivity form connected faces, i.e. similar to the "Laplacian umbrellas"
        con_V2V = con_vertex_vertex_f(F,V)        
    else # e.g. Triangles, which can use simple nodal connectivity
         # Get the unique mesh edges
         E_uni = meshedges(F;unique_only=true) 

         # Get the vertex-to-vertex connectivity, i.e. the "Laplacian umbrellas"
         con_V2V = con_vertex_vertex(E_uni,V)
    end

    NV = vertexnormal(F,V) # The vertex normal directions
    nz = Vec{3,Float64}(0.0,0.0,1.0) # A z-axis vector

    K1 = Vector{Float64}(undef,length(V)) # Allocate first principal curvature
    K2 = Vector{Float64}(undef,length(V)) # Allocate second principal curvature
    U1 = Vector{Vec3{Float64}}(undef,length(V)) # Allocate first principal curvature vector
    U2 = Vector{Vec3{Float64}}(undef,length(V)) # Allocate second principal curvature vector
    for q in eachindex(V)
        n = NV[q] # The current vertex normal
        Q = rotation_between(n,nz) # The rotation between the current normal and the z-axis
        ind = con_V2V[q]        
  
        if growsteps>1            
            for _ in 1:(growsteps-1)
                ind = unique(reduce(vcat,con_V2V[ind]))
                ind = ind[ind.!=q]
            end
        end        
        
        vr = [Q*(v-V[q]) for v in V[ind]] # Rotate point set to a 2D problem
  
        # Set up polynomial fit
        T = Matrix{Float64}(undef,(length(ind),5))
        w = Vector{Float64}(undef,length(ind))
        for i = 1:length(ind)
            T[i,:] = [vr[i][1],vr[i][2],vr[i][1]^2,vr[i][1]*vr[i][2],vr[i][2]^2]
            w[i] = vr[i][3]
        end     

        try
            a = T\w  # x = A\B solves the system of linear equations A*x = B

            E = 1.0 + a[1]^2
            F = a[1]*a[2]
            G = 1.0 + a[2]^2
            d = sqrt(a[1]^2+1.0+a[2]^2)
            e = (2.0*a[3]) / d
            f =       a[4] / d
            g = (2.0*a[5]) / d
    
            S = -[e f; f g]/[E F; F G]
            k,u = eigen(S) # Eigen decomposition to get first/second eigenvalue and vectors
       
            # Store derived quantities
            K1[q] = k[2]
            K2[q] = k[1]
            U1[q] = Q'*Vec3{Float64}(u[1,2],u[2,2],0.0)
            U2[q] = Q'*Vec3{Float64}(u[1,1],u[2,1],0.0)             
        catch
            @warn "SINGULAR"
            K1[q] = NaN
            K2[q] = NaN
            U1[q] = Vec3{Float64}(NaN,NaN,NaN)
            U2[q] = Vec3{Float64}(NaN,NaN,NaN)
        end
    end

    H = 0.5 * (K1.+K2) # Mean curvature
    G = K1.*K2 # Gaussian curvature

    return K1,K2,U1,U2,H,G
end
    
function mesh_curvature_polynomial(M::GeometryBasics.Mesh) 
        return mesh_curvature_polynomial(faces(M),coordinates(M))
end

"""
    separate_vertices(F::Array{NgonFace{N, Int}, 1},V::Array{Point{M, T}, 1}) where N where M where T<:Real
    separate_vertices(M::GeometryBasics.Mesh)

This function takes the input mesh defined by the faces `F` and vertices `V` and
separates any shared vertices. It does this by giving each face its own set of 
unshared vertices. Note that any unused points are not returned in the output 
point array `Vn`. Indices for the mapping are not created here but can simply be
obtained using `reduce(vcat,F)`.
"""
function separate_vertices(F::Union{Vector{NgonFace{N, TF}},Vector{Element{N, TF}}},V::Vector{Point{ND,TV}}; scaleFactor = nothing) where N where TF<:Integer where ND where TV<:Real
    Vn = Vector{eltype(V)}(undef,length(F)*N)
    Fn = Vector{eltype(F)}(undef,length(F))
    for (i,f) in enumerate(F)   
        Fn[i] = (eltype(F))( (i-1)*N .+ (1:N) )       
        if isnothing(scaleFactor)
            Vn[Fn[i]] = V[f]
        else
            vf = mean(V[f])    
            Vn[Fn[i]] = scaleFactor * (V[f] .- vf) .+ vf
        end
    end
    return Fn,Vn
end

function separate_vertices(M::GeometryBasics.Mesh; scaleFactor = nothing)
    F = faces(M)
    V = coordinates(M)
    Fn,Vn = separate_vertices(F,V; scaleFactor = scaleFactor)    
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
    evenly_sample(V::Vector{Point{ND,TV}}, n::Int; rtol = 1e-8, niter = 1) where ND where TV<:Real

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
function evenly_sample(V::Vector{Point{ND,TV}}, n::Int; rtol=1e-8, niter=1, spline_order=4, close_loop=false) where ND where TV<:Real

    S,_,D = make_geospline(V; rtol=rtol, niter=niter, spline_order=spline_order, close_loop=close_loop)

    # Even range for curve distance 
    if close_loop 
        l_end = D - D/n
    else
        l_end = D
    end
    l = range(0.0, l_end, n)     

    return S.(l) # Evaluate interpolator at even distance increments
end

function make_geospline(V::Vector{Point{ND,TV}}; rtol = 1e-8, niter = 10, spline_order=4, close_loop=false) where ND where TV<:Real
    LL = curve_length(V) # Initialise as along curve (multi-linear) distance
    if close_loop
        D = last(LL) + norm(V[1]-V[end])
        bc = BSplineKit.Periodic(D) # Use periodic bc for closed curves
    else
        D = last(LL)  
        bc = BSplineKit.Natural() # Otherwise use natural
    end
    S = BSplineKit.interpolate(LL, deepcopy(V), BSplineOrder(spline_order), bc) # Create interpolator

    L = zeros(eltype(LL),length(LL)) # Initialise spline length vector
    @inbounds for _ in 1:niter
        dS = BSplineKit.Derivative() * S  # spline derivative        
        @inbounds for i in 2:lastindex(LL) 
            # Compute length of segment [i-1, i]   
            L[i] = L[i - 1] + integrate_segment_(dS,LL[i-1], LL[i],rtol)    
        end
        
        if close_loop 
            D = last(L) + integrate_segment_(dS,LL[end], D,rtol)           
            bc = BSplineKit.Periodic(D)                   
        else
            D = last(L)
        end      
        S = BSplineKit.interpolate(L, deepcopy(V), BSplineOrder(spline_order), bc) # Create interpolator
        LL = L
    end
    return S,L,D
end

function integrate_segment_(dS,l1,l2,rtol)
    segment_length, _ = quadgk(l1, l2; rtol) do t
        norm(dS(t))  # integrate |S'(t)| in segment [i, i + 1]
    end    
    return segment_length
end

function evenly_space(V::Vector{Point{ND,TV}}, pointSpacing=nothing; rtol = 1e-8, niter = 1, spline_order=4, close_loop=false, must_points=nothing) where ND where TV<:Real
    if isnothing(pointSpacing)
        pointSpacing = pointspacingmean(V)
    end    

    if !isnothing(must_points)
        sort!(unique!(must_points)) # sort the indices
        if first(must_points)==1 # Check if first is 1
            popfirst!(must_points) # remove start, is already a must point
        end
    end

    if isnothing(must_points) || isempty(must_points)
        n = 1 + ceil(Int,last(curve_length(V; close_loop=close_loop))/pointSpacing)
        Vn = evenly_sample(V,n; rtol=rtol, niter=niter, spline_order=spline_order, close_loop=close_loop)
    else
        # Check must point set
        m = length(V)        
        if !close_loop && last(must_points) != m # Check if last needs to be added 
            push!(must_points,m)
        else
            push!(must_points,1) # Add first to close over curve
        end    

        # Construct length parameterised spline
        S,L,D = make_geospline(V; rtol=rtol, niter=niter, spline_order=spline_order, close_loop=close_loop)

        Vn = Vector{eltype(V)}()
        l1 = 0.0 # Effectively makes first point a must point
        for i in eachindex(must_points)          
            j = must_points[i]          
            if j==1 # i.e. last step when close_loop == true
                l2 = D  
            else
                l2 = L[j]
            end              
            n = 1 + ceil(Int,(l2-l1)/pointSpacing)
            l = range(l1,l2,n)
            Vn_now = S.(l)                
                        
            if i==1 # Append all for first segment
                append!(Vn,Vn_now)  
            elseif j==1 # last step when close_loop == true         
                append!(Vn,Vn_now[2:end-1])
            else # Append 2nd up to last for intermediate segments
                append!(Vn,Vn_now[2:end])  
            end
            l1 = deepcopy(l2)
        end
    end
    return Vn
end

"""
    invert_faces(F::Vector{NgonFace{N, TF}, 1}) where N where TF<:Integer

Flips face orientations.

# Description

This function inverts the faces in `F`, such that the face normal will be 
flipped, by reversing the node order for each face. 
"""
function invert_faces(F::Vector{NgonFace{N, TF}}) where N where TF<:Integer     
    return map(f-> reverse(f),F)     
end

function invert_faces!(F::Vector{NgonFace{N, TF}}) where N where TF<:Integer     
    for (i,f) in enumerate(F)
        F[i] = reverse(f)
    end    
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
    return RotMatrix3{Float64}(V*D*U')
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

    # Determine rotation between sections 
    n = normalizevector(Vc[2]-Vc[1]) # Curve start vector 
    n3_1 = facenormal([collect(1:length(V1))],V1)[1] # Curve normal vector
    if dot(n,n3_1)<0 # Check if anti-aligned
        n3_1 = -n3_1 # Flip vector
        V1 = circshift(reverse(V1),1) # Reverse curve order
        V2 = circshift(reverse(V2),1)
    end

    # Centre start/end sections around means (should be curve start/end points of curve)
    V1b = [v.-Vc[1] for v in V1] 
    V2b = [v.-Vc[end] for v in V2] 
    
    n1_1 = normalizevector(V1b[1]) 
    n2_1 = normalizevector(cross(n3_1,n1_1))
    n1_1 = normalizevector(cross(n2_1,n3_1))
    S1p = mapreduce(permutedims,vcat,[n1_1,n2_1,n3_1])

    n = normalizevector(Vc[end]-Vc[end-1])
    n3_2 = facenormal([collect(1:length(V2))],V2)[1] 
    # if dot(n,n3_2)<0
    #     n3_2 = -n3_2
    # end
    
    n1_2 = normalizevector(V2b[1]) 
    n2_2 = normalizevector(cross(n3_2,n1_2))
    n1_2 = normalizevector(cross(n2_2,n3_2))        
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
        else # Rotate and position intermediate section 
            if q == nc
                n3 = normalizevector(Vc[q]-Vc[q-1])                        
            else
                n3 = normalizevector(normalizevector(Vc[q]-Vc[q-1]) .+ normalizevector(Vc[q+1]-Vc[q]))                        
            end
            n2 = normalizevector(cross(n3,n1p))
            n1 = normalizevector(cross(n2,n3))
            S2 = mapreduce(permutedims,vcat,[n1,n2,n3])

            Q12 = RotMatrix3{Float64}(S2\S1p)
            V[ind] = [(Q12*v).+Vc[q] for v in V[ind]]
            n1p = n1
        end   

        if q == nc # Once here, a potential rotational mismatch needs to be resolved
            V[ind] = V2 # Overwrite end with desired end

            # Determine rotation between last and actual end section
            Q_fix = RotMatrix3{Float64}(S2\S2p) 
            
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
    revolvecurve(Vc::Vector{Point{ND,TV}}; extent = 2.0*pi, direction=:positive, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing,close_loop=true,face_type=:quad)  where ND where TV<:Real   

Revolves curves to build surfaces 

# Description

This function rotates the curve `Vc` by the angle `extent`, in the direction 
defined by `direction` (`:positive`, `:negative`, `:both`), around the vector 
`n`, to build the output mesh defined by the faces `F` and vertices `V`. 
"""
function revolvecurve(Vc::Vector{Point{ND,TV}}; extent = 2.0*pi, direction=:positive, n=Vec{3, Float64}(0.0,0.0,1.0),num_steps=nothing, periodicity=(false,false),face_type=:quad)  where ND where TV<:Real   
    
    # Compute num_steps from curve point spacing
    if isnothing(num_steps)
        rMax = 0.0
        for v in Vc
            rNow = dot(normalizevector(cross(cross(n,v),n)),v)
            if !isnan(rNow)
                rMax = max(rMax,rNow)
            end
        end
        num_steps = ceil(Int,(rMax*extent)/pointspacingmean(Vc))        
    end
    
    # Set up angle range
    if direction == :positive # Positive direction
        θ_range = range(0,extent,num_steps)               
    elseif direction == :negative # Negative direction
        θ_range = range(-extent,0,num_steps) # Positive range         
    elseif direction == :both # Both positive and negative directions 
        θ_range = range(-extent/2,extent/2,num_steps) # Positive range   
    else
        throw(ArgumentError("$direction is not a valid direction, Use :positive, :in, or :both.")) 
    end

    V = Vector{eltype(Vc)}()
    for θ in θ_range
        Q = AngleAxis(θ, n[1], n[2], n[3])
        Vn = [Q*v for v in Vc]
        append!(V,Vn)
    end    
   
    F = grid2surf(V,num_steps; face_type=face_type, periodicity=periodicity)

    return F,V
end


"""
    batman(n::Int)
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
function batman(n::Int; symmetric = false, dir=:acw)
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

    if symmetric    
        if iseven(n)
            m = round(Int,n/2+1)
        else
            m = ceil(Int,n/2)+1 
        end
        V = evenly_sample([Point{3, Float64}(x[i]/22,y[i]/22,0.0) for i in eachindex(x)],m)
        V2 = [Point{3, Float64}(-V[i][1],V[i][2],0.0) for i in length(V)-1:-1:2]
        append!(V,V2)    
    else
        V = [Point{3, Float64}(x[i]/22,y[i]/22,0.0) for i in eachindex(x)]    
        V2 = [Point{3, Float64}(-V[i][1],V[i][2],0.0) for i in length(V)-1:-1:2]
        append!(V,V2) 
        V = evenly_sample(V,n)
    end
    if dir==:cw
        reverse!(V)   
        circshift!(V,1)     
    end
    V = V.-mean(V,dims=1)
    return V
end

"""
    tridisc(r=1.0,n=0; ngon=6, method = :linear, orientation=:up)

# Description

Generates the faces `F` and vertices `V` for a triangulated disc (circle). The 
algorithm starts with a triangulated hexagon (obtained if `n=0`) and uses 
iterative subtriangulation, and uses iterative subdivision (and pushing of 
boundary points to circular boundary) to obtain the final mesh. The subdivision
`method` is an optional input, and is either `:Loop` (default) or `:linear`. 
Lastly the optional input `orientation`, which can be `:up` or `:down` sets the 
face normal direction. 
"""
function tridisc(r=1.0,n=0; ngon=6, method=:Loop, orientation=:up)
    if  !in(orientation,(:up,:down))        
        throw(ArgumentError("Orientation not supported. Use :up or :down"))
    end 

    # Create a triangulated hexagon
    V = circlepoints(r,ngon)
    push!(V,Point{3,Float64}(0.0,0.0,0.0))
    if orientation==:up
        F = [TriangleFace{Int}(i,mod1(i+1,ngon),ngon+1) for i in 1:ngon]
    else#if orientation==:down
        F = [TriangleFace{Int}(i,ngon+1,mod1(i+1,ngon)) for i in 1:ngon]    
    end

    # Refine mesh n times
    if n>0
        @inbounds for _ = 1:n            
            nv = length(V)             
            F,V = subtri(F,V,1; method = method)
            indBoundary = unique(reduce(vcat,boundaryedges(F)))
            @inbounds for i in indBoundary
                if i>nv || method == :Loop
                    V[i] *= r/norm(V[i])
                end
            end            
        end        
    end
    return F,V
end

"""
    triplate(xSpan,ySpan,pointSpacing::T) where T <: Real    

# Description
Generates a triangulated mesh for a plate.
"""
function triplate(plateDim,pointSpacing::T; orientation=:up) where T <: Real 
    if  !in(orientation,(:up,:down))        
        throw(ArgumentError("Orientation not supported. Use :up or :down"))
    end 

    F, V = gridpoints_equilateral((-plateDim[1]/2,plateDim[1]/2),(-plateDim[2]/2,plateDim[2]/2),pointSpacing; return_faces = true, rectangular=true, force_equilateral=false)
    if orientation == :down
        return invert_faces(F), V
    else
        return F,V
    end
end

"""
    quaddisc(r,n; method = :Catmull_Clark, orientation=:up)

# Description

Generates the faces `F` and vertices `V` for a quadrangulated (circle). The 
algorithm starts with 12 quadrilateral faces and an octagon boundary  
(obtained if `n=0`), and uses iterative subdivision (and pushing of boundary 
points to circular boundary) to obtain the final mesh. The subdivision
`method` is an optional input, and is either `:Catmull_Clark` (default) or 
`:linear`. Lastly the optional input `orientation`, which can be `:up` or 
`:down` sets the face normal direction. 
"""
function quaddisc(r,n; method = :Catmull_Clark, orientation=:up)

    if  !in(orientation,(:up,:down))        
        throw(ArgumentError("Orientation not supported. Use :up or :down"))
    end 

    # Create disc mesh 
    V = circlepoints(r,8)
    V = append!(V, circlepoints(r/2,8)) 
    V = push!(V, Point{3,Float64}(0.0,0.0,0.0))
    F = Vector{QuadFace{Int}}(undef,12)

    # Add outer ring
    for i = 1:8
        F[i] = QuadFace{Int}(i, mod1(i+1,8), mod1(i+1+8,8)+8, mod1(i+8,16))
    end

    # Add core
    for i = 1:4
        j = 9 + (i-1)*2
        F[i+8] = QuadFace{Int}(j, mod1(j+1+8,8)+8, mod1(j+2+8,8)+8, 17)
    end
    
    # Invert if needed
    if orientation==:down
        F = invert_faces(F)    
    end
    
    # Refine
    for _ = 1:1:n
        nv = length(V)             
        F,V = subquad(F,V,1; method=method) 
        indBoundary = unique(reduce(vcat,boundaryedges(F)))
        @inbounds for i in indBoundary
            if i>nv || method == :Catmull_Clark
                V[i] *= r/norm(V[i])
            end
        end        
    end
    return F,V
end


"""
    regiontrimesh(VT,R,P)

# Description

Generates a multi-region triangle mesh for the input regions. The boundary 
curves for all regions are contained in the tuple `VT`. Each region to be meshed
is next defined using a tuple `R` containing indices into the curve typle `VT`. 
If an entry in `R` contains only one index then the entire curve domain is 
meshed. If `R` contains multiple indices then the first index is assumed to be 
for the outer boundary curve, while all subsequent indices are for boundaries 
defining holes in this region. 
"""
function regiontrimesh(VT,R,P)
    V = Vector{Point{3,Float64}}() # eltype(VT)()
    F = Vector{TriangleFace{Int}}()
    C = Vector{Float64}()
    
    for q in eachindex(R)        
        r = R[q] # The curve indices for the current region       
        pointSpacing = P[q] # Region point spacing
        Vn = reduce(vcat,VT[r]) # The current curve point set 
        
        # Creating constraints
        constrained_segments = Vector{Vector{Vector{Int}}}()
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
        zMean = mean([v[3] for v in Vn])
        Vn = append!(Vn,Vg)        
        Vn = [Point{3,Float64}(v[1],v[2],zMean) for v in Vn] # Force zero z-coordinate
        
        # Check for unique points 
        np = length(Vn)
        Vn,_,indMap = mergevertices(Vn; pointSpacing=pointSpacing)
        if length(Vn)<np
            for (i,c) in enumerate(constrained_segments)
                constrained_segments[i] = [indMap[cc] for cc in c]
            end
        end

        # Initial triangulation 
        constrained_segments_ori = deepcopy(constrained_segments) # Clone since triangulate can add new constraint points
        TRn = triangulate(Vn; boundary_nodes=constrained_segments, delete_ghosts=true)
        Fn = [TriangleFace{Int}(tr) for tr in each_solid_triangle(TRn)] 
        Vn = get_points(TRn)
        # Fn,Vn,indFix = remove_unused_vertices(Fn,Vn)
        # constrained_segments = [[indFix[c[1]]] for c in constrained_segments_ori] # Fix indices after point removal 
    
        # Check if new boundary points were introduced and remove if needed 
        Eb = boundaryedges(Fn)
        indB = unique(reduce(vcat,Eb))
        indConstrained = reduce(vcat,reduce(vcat,constrained_segments_ori))
        indRemove = setdiff(indB,indConstrained)
        if !isempty(indRemove)    
            Vn,indFix = removepoints(Vn,indRemove)
            constrained_segments = [[indFix[c[1]]] for c in constrained_segments_ori]
    
            # Redo triangulation after points have been removed
            constrained_segments_ori = deepcopy(constrained_segments) # Clone since triangulate can add new constraint points
            TRn = triangulate(Vn; boundary_nodes=constrained_segments,delete_ghosts=true)
            Fn = [TriangleFace{Int}(tr) for tr in each_solid_triangle(TRn)] 
            Vn = get_points(TRn)
            constrained_segments = constrained_segments_ori
        end
    
        # Remove 3 and 4 connected points
        E_uni = meshedges(Fn; unique_only=false)
        con_V2V = con_vertex_vertex(E_uni,Vn)
        nCon = map(length,con_V2V)        
        indLowCon = findall(nCon.>0 .&& nCon.<5)
        indConstrained = reduce(vcat,reduce(vcat,constrained_segments))
        indRemove = setdiff(indLowCon,indConstrained)     
        if !isempty(indRemove)
            Vn,indFix = removepoints(Vn,indRemove)
            constrained_segments = [[indFix[c[1]]] for c in constrained_segments]
    
            # Redo triangulation after points have been removed
            constrained_segments_ori = deepcopy(constrained_segments) # Clone since triangulate can add new constraint points
            TRn = triangulate(Vn; boundary_nodes=constrained_segments,delete_ghosts=true)
            Fn = [TriangleFace{Int}(tr) for tr in each_solid_triangle(TRn)] 
            Vn = get_points(TRn)
            constrained_segments = constrained_segments_ori
        end    
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
    V,_,indMap = mergevertices(V; pointSpacing = mean(P))
    indexmap!(F,indMap)    
    
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
    Fs,Vs = separate_vertices(F,V) # separate_vertices
    for (i,f) in enumerate(Fs) 
        vm = mean(Vs[f],dims=1)
        if length(s)==length(F)
            Vs[f] = s[i].*(Vs[f].-vm).+vm    
        elseif length(s)==length(V)
            Vs[f] = s[F[i]].*(Vs[f].-vm).+vm    
        elseif isone(length(s))
            Vs[f] = s.*(Vs[f].-vm).+vm    
        end    
    end
    return Fs,Vs
end

"""
    subcurve(V,n)

Adds `n` points between each curve point.  

# Description

This function adds `n` points between each current point on the curve `V`.
"""
function subcurve(V::Vector{Point{ND,TV}}, n; close_loop=false) where ND where TV<:Real    
    m = length(V)
    
    if close_loop == true
        Vn = Vector{Point{ND,TV}}(undef,m+m*n)    
    else
        Vn = Vector{Point{ND,TV}}(undef,m+(m-1)*n)    
    end
    
    @inbounds for q = 1:lastindex(V)-1                
        vn = range(V[q],V[q+1],n+2)
        i1 = 1+(q-1)*(n+1)
        Vn[i1:i1+n] = vn[1:end-1]
    end

    if close_loop == true
        vn = range(V[end],V[1],n+2)
        i1 = 1+(m-1)*(n+1)
        Vn[i1:i1+n] = vn[1:end-1]
    else
        Vn[end] = V[end]
    end
    return Vn
end


"""
    dualclad(F::Vector{NgonFace{N, TF}},V::Vector{Point{ND,TV}},s::T; connectivity=:face) where N where TF<:Integer where ND where TV<:Real where T<:Real

Returns a surface conforming dual lattice

# Description
 
"""
function dualclad(F::Vector{NgonFace{N, TF}},V::Vector{Point{ND,TV}},s; connectivity=:face) where N where TF<:Integer where ND where TV<:Real
    Fs,Vs = scalesimplex(F,V,s) # Scaled faces

    # All non-unique mesh edges 
    E = meshedges(F;unique_only=false) 
    E_Fs = meshedges(Fs;unique_only=false)

    if connectivity == :face 
        # Compute dict to find partner edges, unique edges as well as reverse indices for unique mapping
        E_sort = sort.(E) 
        d = Dict{eltype(E),Vector{Int}}() # Use dict to keep track of used values    
        Eu = Vector{eltype(E)}()
        indReverse = Vector{Int}(undef,length(E)) 
        j=0
        for (i,e) in enumerate(E_sort)          
            if !haskey(d, e)    
                j+=1 # Increment counter        
                d[e]= [i]  
                indReverse[i] = j  # Store inverse index    
                push!(Eu, E[i]) #Grow unique set      
            else
                push!(d[e],i)    
                indReverse[i] = indReverse[d[e][1]]      
            end
        end  

        # Check for boundary faces 
        count_E2F = count_edge_face(F,Eu,indReverse)
        indBoundary = findall(isone.(count_E2F))
        Eb = Eu[indBoundary]

        indX = zeros(Int,length(Eu))
        indX[indBoundary]=1:length(indBoundary)
        indX_E = indX[indReverse]

        if !isempty(Eb)
            Ebs,Vbs = scalesimplex(Eb,V,s)
            Ebs,Vbs,_ = remove_unused_vertices(Ebs,Vbs)
            Ebs = [e.+length(Vs) for e in Ebs]
            append!(Vs,Vbs) # Append boundary edge points 
        end

        Fq = Vector{QuadFace{Int}}(undef,length(Eu))
        for (i,e) in enumerate(Eu)
            ii = d[sort(e)]
            if length(ii)==2
                Fq[i] = (E_Fs[ii[1]][2],E_Fs[ii[1]][1],E_Fs[ii[2]][2],E_Fs[ii[2]][1])
            else
                ii_b = indX_E[ii[1]]
                Fq[i] = (E_Fs[ii[1]][2],E_Fs[ii[1]][1],Ebs[ii_b][1],Ebs[ii_b][2])
            end
        end
    elseif connectivity == :edge
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

        Es,V_Es = scalesimplex(Eu,V,s)
        Es = [e.+length(Vs) for e in Es]
        append!(Vs,V_Es) # Append boundary edge points 
        
        Fq = Vector{QuadFace{Int}}(undef,length(E_Fs))
        for i in eachindex(E_Fs)
            ii = indReverse[i]
            if E[i][1] == Eu[ii][1] # A first occurance edge, order matches
                Fq[i] = (E_Fs[i][2],E_Fs[i][1],Es[ii][1],Es[ii][2])
            else # A second occurance edge, needs inversion
                Fq[i] = (E_Fs[i][2],E_Fs[i][1],Es[ii][2],Es[ii][1])
            end
        end
    end

    return Fs,Fq,Vs
end


"""
    tet2hex(E::Vector{Tet4{T}},V::Vector{Point{ND,TV}}) where T<:Integer where ND where TV<:Real

Converts tetrahedra to hexahedra

# Description
This function converts the input tetrahedra defined by the element set `E` and the 
vertex set `V` to a set of hexahedral elements `Eh` with vertices `Vh`. The conversion
involves a splitting of each tetrahedron into 4 hexahedra. 
"""
function tet2hex(E::Vector{Tet4{T}},V::Vector{Point{ND,TV}}) where T<:Integer where ND where TV<:Real
    # Non-unique tet element faces
    Ft = element2faces(E)      
    
    # Unique faces and inverse mapping 
    Ftu,_,inv_Ft = gunique(Ft; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

    # Non-unique structured (element by element) tet element edges
    Et = Vector{LineFace{T}}(undef,length(E)*6)
    for (i,e) in enumerate(E)         
        ii = 1 + (i-1)*6
        Et[ii  ] = LineFace{T}(e[1],e[2])
        Et[ii+1] = LineFace{T}(e[2],e[3])
        Et[ii+2] = LineFace{T}(e[3],e[1])
        Et[ii+3] = LineFace{T}(e[1],e[4])
        Et[ii+4] = LineFace{T}(e[2],e[4])
        Et[ii+5] = LineFace{T}(e[3],e[4])
    end

    #Unique edges and inverse mapping 
    Etu,_,inv_Et = gunique(Et; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
    
    # Create new coordinates
    Vc = simplexcenter(E,V) # Element centres
    Vf = simplexcenter(Ftu,V) # Face centres
    Ve = simplexcenter(Etu,V) # Edge centres   
    Vh = [V;Vc;Vf;Ve] # Append vertices

    # Create the 4 corner hexahedra for each tetrahedron
    offset_Vc = length(V)
    offset_F = offset_Vc + length(Vc)
    offset_E = offset_F + length(Vf)
    inv_Ft = inv_Ft .+ offset_F
    inv_Et = inv_Et .+ offset_E
    Eh = Vector{Hex8{T}}(undef,length(E)*4) # Allocate hexahedral element vector, one hex for each of the 4 nodes per element  
    for i = eachindex(E)
        i_e = 1 + (i-1)*6 # index of first edge
        i_f = 1 + (i-1)*4 # index of first face (which happens to also be the first hex element for tetrahedra)        
        i_vc = i+offset_Vc
        Eh[i_f  ] = Hex8{T}(      E[i][1], inv_Et[i_e  ], inv_Ft[i_f  ], inv_Et[i_e+2],
                            inv_Et[i_e+3], inv_Ft[i_f+1],          i_vc, inv_Ft[i_f+3])
        Eh[i_f+1] = Hex8{T}(      E[i][2], inv_Et[i_e+1], inv_Ft[i_f  ], inv_Et[i_e  ],
                            inv_Et[i_e+4], inv_Ft[i_f+2],          i_vc, inv_Ft[i_f+1])
        Eh[i_f+2] = Hex8{T}(      E[i][3], inv_Et[i_e+2], inv_Ft[i_f  ], inv_Et[i_e+1],
                            inv_Et[i_e+5], inv_Ft[i_f+3],          i_vc, inv_Ft[i_f+2])                             
        Eh[i_f+3] = Hex8{T}(      E[i][4], inv_Et[i_e+4], inv_Ft[i_f+1], inv_Et[i_e+3],
                            inv_Et[i_e+5], inv_Ft[i_f+2],          i_vc, inv_Ft[i_f+3])                             
    end    
    return Eh,Vh
end

"""
    element2faces(E::Vector{Element{N,T}}) where N where T 

Returns element faces

# Description
This function computes the faces for the input elements defined by `E`. The elements
should be Vectors consisting of `Tet4`, `Hex8` elements. 
"""
function element2faces(E::Vector{Element{N,T}}) where N where T 
    element_type = eltype(E)
    if element_type <: Tet4{T}
        nf = 4
        F = Vector{TriangleFace{T}}(undef,length(E)*nf)
        for (i,e) in enumerate(E)
            ii = 1 + (i-1)*nf
            F[ii  ] = TriangleFace{T}(e[3],e[2],e[1])
            F[ii+1] = TriangleFace{T}(e[1],e[2],e[4])
            F[ii+2] = TriangleFace{T}(e[2],e[3],e[4])
            F[ii+3] = TriangleFace{T}(e[3],e[1],e[4])
        end
    elseif element_type <: Tet10{T}
        nf = 4        
        F = Vector{NgonFace{6,Int}}(undef,length(E)*nf)
        for (i,e) in enumerate(E)
            ii = 1 + (i-1)*nf 
            F[ii  ] = NgonFace{6,Int}(e[3],e[6],e[2],e[5 ],e[1],e[7 ])
            F[ii+1] = NgonFace{6,Int}(e[1],e[5],e[2],e[9 ],e[4],e[8 ])
            F[ii+2] = NgonFace{6,Int}(e[2],e[6],e[3],e[10],e[4],e[9 ])
            F[ii+3] = NgonFace{6,Int}(e[3],e[7],e[1],e[8 ],e[4],e[10])
        end
    elseif element_type <: Hex8{T}
        nf = 6
        F = Vector{QuadFace{T}}(undef,length(E)*nf)
        for (i,e) in enumerate(E)            
            ii = 1 + (i-1)*nf
            F[ii  ] = QuadFace{T}(e[4],e[3],e[2],e[1]) # Top
            F[ii+1] = QuadFace{T}(e[5],e[6],e[7],e[8]) # Bottom
            F[ii+2] = QuadFace{T}(e[1],e[2],e[6],e[5]) # Side 1
            F[ii+3] = QuadFace{T}(e[8],e[7],e[3],e[4]) # Side 2
            F[ii+4] = QuadFace{T}(e[2],e[3],e[7],e[6]) # Front
            F[ii+5] = QuadFace{T}(e[5],e[8],e[4],e[1]) # Back
        end
    elseif element_type <: Penta6{T}
        nfq = 3
        Fq = Vector{QuadFace{T}}(undef,length(E)*nfq)
        for (i,e) in enumerate(E)
            ii = 1 + (i-1)*nfq
            Fq[ii  ] = QuadFace{T}(e[1],e[2],e[5],e[4]) # Side 1
            Fq[ii+1] = QuadFace{T}(e[2],e[3],e[6],e[5]) # Side 2
            Fq[ii+2] = QuadFace{T}(e[3],e[1],e[4],e[6]) # Side 3            
        end

        nft = 2
        Ft = Vector{TriangleFace{T}}(undef,length(E)*nft)
        for (i,e) in enumerate(E)
            ii = 1 + (i-1)*nft
            Ft[ii  ] = TriangleFace{T}(e[3],e[2],e[1]) # Bottom
            Ft[ii+1] = TriangleFace{T}(e[4],e[5],e[6]) # Top
        end
        F = (Ft,Fq) # Collect faces in tuple
    elseif element_type <: Rhombicdodeca14{T}
        F0 = Vector{QuadFace{Int}}(undef,12)
        F0[ 1] = QuadFace{Int}( 1, 10, 5,  9)
        F0[ 2] = QuadFace{Int}( 2, 11, 6, 10)
        F0[ 3] = QuadFace{Int}( 3, 12, 7, 11)
        F0[ 4] = QuadFace{Int}( 4,  9, 8, 12)
        F0[ 5] = QuadFace{Int}( 5, 10, 6, 14)
        F0[ 6] = QuadFace{Int}( 6, 11, 7, 14)
        F0[ 7] = QuadFace{Int}( 7, 12, 8, 14)
        F0[ 8] = QuadFace{Int}( 8,  9, 5, 14)
        F0[ 9] = QuadFace{Int}( 1, 13, 2, 10)
        F0[10] = QuadFace{Int}( 2, 13, 3, 11)
        F0[11] = QuadFace{Int}( 3, 13, 4, 12)
        F0[12] = QuadFace{Int}( 4, 13, 1,  9)

        F = Vector{QuadFace{Int}}(undef,length(E)*length(F0))
        for (i,e) in enumerate(E) 
            j = 1+(i-1)*length(F0)
            F[j:j+length(F0)-1] = [QuadFace{Int}(view(e,f)) for f in F0]    
        end
    elseif element_type <: Truncatedocta24{T}
        # Hexagonal faces
        F01 = Vector{NgonFace{6,Int}}(undef,8)
        F01[ 1] = NgonFace{6,Int}( 6, 18,  1, 13,  3, 15)
        F01[ 2] = NgonFace{6,Int}( 9, 21,  5, 17, 18,  6)
        F01[ 3] = NgonFace{6,Int}(11, 23,  8, 20, 21,  9)
        F01[ 4] = NgonFace{6,Int}(15,  3,  2, 14, 23, 11)
        F01[ 5] = NgonFace{6,Int}(13,  1,  7, 19,  4, 16)
        F01[ 6] = NgonFace{6,Int}(17,  5, 10, 22, 19,  7)
        F01[ 7] = NgonFace{6,Int}(20,  8, 12, 24, 22, 10)
        F01[ 8] = NgonFace{6,Int}(14,  2, 16,  4, 24, 12)

        # Quadrilateral faces
        F02 = Vector{QuadFace{Int}}(undef,6)
        F02[ 1] = QuadFace{Int}(2,  3, 13, 16)
        F02[ 2] = QuadFace{Int}(17,  7,  1, 18)
        F02[ 3] = QuadFace{Int}(20, 10,  5, 21)
        F02[ 4] = QuadFace{Int}(14, 12,  8, 23)
        F02[ 5] = QuadFace{Int}(11,  9,  6, 15)
        F02[ 6] = QuadFace{Int}(4, 19, 22, 24)

        F1 = Vector{NgonFace{6,Int}}(undef,length(E)*length(F01))
        for (i,e) in enumerate(E) 
            j = 1+(i-1)*length(F01)
            F1[j:j+length(F01)-1] = [NgonFace{6,Int}(view(e,f)) for f in F01]    
        end

        F2 = Vector{QuadFace{Int}}(undef,length(E)*length(F02))
        for (i,e) in enumerate(E) 
            j = 1+(i-1)*length(F02)
            F2[j:j+length(F02)-1] = [QuadFace{Int}(view(e,f)) for f in F02]    
        end
        F = (F1,F2) # Collect faces in tuple
    else
        throw(ArgumentError("$element_type not supported. Supported types are Hex8, Tet4, and Penta6."))
    end
    return F
end

"""
    subhex(E::Vector{Hex8{T}},V::Vector{Point{ND,TV}},n::Int; direction=0) where T<:Integer where ND where TV<:Real

Split hexahedral elements

# Description
This function splits the hexahedral elements defined by the elements `E` and 
vertices `V`. Splitting is done `n` times as requested. By default the splitting 
occurs in all direction (corresponding to the default `direction=0`). If instead
`direction` is set to 1, 2, or 3, then the splitting only occur in the first, second
or third local element direction respectively. Note that this direction depends 
on node order used. For a hexahedron where by nodes 1:4 are for the bottom, and 
nodes 5:8 are for the top of the element then the directions 1, 2, and 3 correspond 
to the x-, y-, and z-direction respectively.  
"""
function subhex(E::Vector{Hex8{T}},V::Vector{Point{ND,TV}},n::Int; direction=0) where T<:Integer where ND where TV<:Real
    if iszero(n)
        return E,V
    elseif isone(n)

        # Non-unique hex element faces
        Ft = element2faces(E)      
        
        # Non-unique structured (element by element) hexbox element edges
        if iszero(direction)
            Et = Vector{LineFace{T}}(undef,length(E)*12)
            for (i,e) in enumerate(E)         
                ii = 1 + (i-1)*12
                Et[ii   ] = LineFace{T}(e[1],e[2])
                Et[ii+1 ] = LineFace{T}(e[2],e[3])
                Et[ii+2 ] = LineFace{T}(e[3],e[4])
                Et[ii+3 ] = LineFace{T}(e[4],e[1])
                Et[ii+4 ] = LineFace{T}(e[5],e[6])
                Et[ii+5 ] = LineFace{T}(e[6],e[7])
                Et[ii+6 ] = LineFace{T}(e[7],e[8])
                Et[ii+7 ] = LineFace{T}(e[8],e[5])
                Et[ii+8 ] = LineFace{T}(e[1],e[5])
                Et[ii+9 ] = LineFace{T}(e[2],e[6])
                Et[ii+10] = LineFace{T}(e[3],e[7])
                Et[ii+11] = LineFace{T}(e[4],e[8])
            end
        elseif isone(direction) # Split in 1st-direction
            Et = Vector{LineFace{T}}(undef,length(E)*4)
            for (i,e) in enumerate(E)         
                ii = 1 + (i-1)*4
                Et[ii   ] = LineFace{T}(e[1],e[2])
                Et[ii+1 ] = LineFace{T}(e[3],e[4])
                Et[ii+2 ] = LineFace{T}(e[5],e[6])
                Et[ii+3 ] = LineFace{T}(e[7],e[8])
            end
        elseif direction==2 # Split in 2nd-direction
            Et = Vector{LineFace{T}}(undef,length(E)*4)
            for (i,e) in enumerate(E)         
                ii = 1 + (i-1)*4
                Et[ii   ] = LineFace{T}(e[2],e[3])
                Et[ii+1 ] = LineFace{T}(e[4],e[1])
                Et[ii+2 ] = LineFace{T}(e[6],e[7])
                Et[ii+3 ] = LineFace{T}(e[8],e[5])
            end
        elseif direction==3 # Split in 3rd-direction
            Et = Vector{LineFace{T}}(undef,length(E)*4)
            for (i,e) in enumerate(E)         
                ii = 1 + (i-1)*4
                Et[ii   ] = LineFace{T}(e[1],e[5])
                Et[ii+1 ] = LineFace{T}(e[2],e[6])
                Et[ii+2 ] = LineFace{T}(e[3],e[7])
                Et[ii+3 ] = LineFace{T}(e[4],e[8])
            end
        end
        #Unique edges and inverse mapping 
        Etu,_,inv_Et = gunique(Et; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        # Create mid-edge coordinates
        Ve = simplexcenter(Etu,V) # Edge centres           

        if iszero(direction)
            # Unique faces and inverse mapping 
            Ftu,_,inv_Ft = gunique(Ft; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
            Vf = simplexcenter(Ftu,V) # Face centres
            Vc = simplexcenter(E,V) # Element centres
            Vh = [V;Vc;Vf;Ve] # Append vertices

            # Create the 8 corner hexahedra for each hexahedron
            offset_Vc = length(V)
            offset_F = offset_Vc + length(Vc)
            offset_E = offset_F + length(Vf)
            inv_Ft = inv_Ft .+ offset_F
            inv_Et = inv_Et .+ offset_E  
            Eh = Vector{Hex8{T}}(undef,length(E)*8) # Allocate hexahedral element vector, one hex for each of the 4 nodes per element  
            for i = eachindex(E)
                i_e = 1 + (i-1)*12 # index of first edge
                i_f = 1 + (i-1)*6 # index of first face 
                i_vc = i + offset_Vc        
                ii = 1 + (i-1)*8 
                Eh[ii  ] = Hex8{T}(        E[i][1], inv_Et[i_e  ], inv_Ft[i_f  ], inv_Et[i_e+3],
                                    inv_Et[i_e+8 ], inv_Ft[i_f+2],          i_vc, inv_Ft[i_f+5] )   

                Eh[ii+1] = Hex8{T}(        E[i][2], inv_Et[i_e+1], inv_Ft[i_f  ], inv_Et[i_e],
                                    inv_Et[i_e+9 ], inv_Ft[i_f+4],          i_vc, inv_Ft[i_f+2] )   

                Eh[ii+2] = Hex8{T}(        E[i][3], inv_Et[i_e+2], inv_Ft[i_f  ], inv_Et[i_e+1],
                                    inv_Et[i_e+10], inv_Ft[i_f+3],          i_vc, inv_Ft[i_f+4] )   

                Eh[ii+3] = Hex8{T}(        E[i][4], inv_Et[i_e+3], inv_Ft[i_f  ], inv_Et[i_e+2],
                                    inv_Et[i_e+11], inv_Ft[i_f+5],          i_vc, inv_Ft[i_f+3] )   
            
                Eh[ii+4] = Hex8{T}(        E[i][5], inv_Et[i_e+7], inv_Ft[i_f+1], inv_Et[i_e+4],
                                    inv_Et[i_e+8 ], inv_Ft[i_f+5],          i_vc, inv_Ft[i_f+2] )   

                Eh[ii+5] = Hex8{T}(        E[i][6], inv_Et[i_e+4], inv_Ft[i_f+1], inv_Et[i_e+5],
                                    inv_Et[i_e+9 ], inv_Ft[i_f+2],          i_vc, inv_Ft[i_f+4] )   

                Eh[ii+6] = Hex8{T}(        E[i][7], inv_Et[i_e+5], inv_Ft[i_f+1], inv_Et[i_e+6],
                                    inv_Et[i_e+10], inv_Ft[i_f+4],          i_vc, inv_Ft[i_f+3] )   
                                    
                Eh[ii+7] = Hex8{T}(        E[i][8], inv_Et[i_e+6], inv_Ft[i_f+1], inv_Et[i_e+7],
                                    inv_Et[i_e+11], inv_Ft[i_f+3],          i_vc, inv_Ft[i_f+5] )  
            end
        elseif isone(direction) # Split in 1st-direction
            Vh = [V;Ve] # Append vertices
            inv_Et = inv_Et .+ length(V)  
            Eh = Vector{Hex8{T}}(undef,length(E)*2) # Allocate hexahedral element vector, one hex for each of the 4 nodes per element  
            for i = eachindex(E)
                i_e = 1 + (i-1)*4 # index of first edge
                i_f = 1 + (i-1)*6 # index of first face 
                ii = 1 + (i-1)*2
                Eh[ii  ] = Hex8{T}( Ft[i_f+4][4],  Ft[i_f+4][3],  Ft[i_f+4][2],  Ft[i_f+4][1], 
                                   inv_Et[i_e+2], inv_Et[i_e+3], inv_Et[i_e+1], inv_Et[i_e+0] )   

                Eh[ii+1] = Hex8{T}( Ft[i_f+5][4],  Ft[i_f+5][3],  Ft[i_f+5][2],  Ft[i_f+5][1], 
                                   inv_Et[i_e+0], inv_Et[i_e+1], inv_Et[i_e+3], inv_Et[i_e+2] )      

            end  
        elseif direction == 2 # Split in 2nd-direction
            Vh = [V;Ve] # Append vertices
            inv_Et = inv_Et .+ length(V)  
            Eh = Vector{Hex8{T}}(undef,length(E)*2) # Allocate hexahedral element vector, one hex for each of the 4 nodes per element  
            for i = eachindex(E)
                i_e = 1 + (i-1)*4 # index of first edge
                i_f = 1 + (i-1)*6 # index of first face 
                ii = 1 + (i-1)*2
                Eh[ii  ] = Hex8{T}( Ft[i_f+2][4],  Ft[i_f+2][3],  Ft[i_f+2][2],  Ft[i_f+2][1], 
                                   inv_Et[i_e+3], inv_Et[i_e+2], inv_Et[i_e+0], inv_Et[i_e+1] )   

                Eh[ii+1] = Hex8{T}( Ft[i_f+3][4],  Ft[i_f+3][3],  Ft[i_f+3][2],  Ft[i_f+3][1],
                                   inv_Et[i_e+1], inv_Et[i_e+0], inv_Et[i_e+2], inv_Et[i_e+3] )   
            end    
        elseif direction == 3 # Split in 3rd-direction
            Vh = [V;Ve] # Append vertices
            inv_Et = inv_Et .+ length(V)  
            Eh = Vector{Hex8{T}}(undef,length(E)*2) # Allocate hexahedral element vector, one hex for each of the 4 nodes per element  
            for i = eachindex(E)
                i_e = 1 + (i-1)*4 # index of first edge
                i_f = 1 + (i-1)*6 # index of first face 
                ii = 1 + (i-1)*2
                Eh[ii  ] = Hex8{T}(   Ft[i_f][4],    Ft[i_f][3],     Ft[i_f][2],     Ft[i_f][1], 
                                   inv_Et[i_e+0], inv_Et[i_e+1], inv_Et[i_e+2], inv_Et[i_e+3] )   

                Eh[ii+1] = Hex8{T}(  Ft[i_f+1][4],   Ft[i_f+1][3],  Ft[i_f+1][2],  Ft[i_f+1][1],
                                   inv_Et[i_e+3], inv_Et[i_e+2], inv_Et[i_e+1], inv_Et[i_e+0] )   

            end   
        end 
        return Eh,Vh
    elseif n>1
        for _ = 1:n
            E,V = subhex(E,V,1; direction=direction)
        end
        return E,V
    else
        throw(ArgumentError("n should be larger than or equal to 0"))
    end
end

"""
    rhombicdodecahedron(r = 1.0)

Creates mesh for rhombicdodecahedron

# Description
This function creates the faces `F` and vertices `V` for a rhombicdodecahedron. 
The size of the shape is set by the width `w` in the xy-plane. 
"""
function rhombicdodecahedron(w = 1.0)

    a = sqrt(2)/2
    b = sqrt(2)/4
    V = Vector{Point{3,Float64}}(undef,14)
    V[ 1] = Point{3,Float64}(-0.0, -0.5,  -b)
    V[ 2] = Point{3,Float64}( 0.5, -0.0,  -b)
    V[ 3] = Point{3,Float64}( 0.0,  0.5,  -b)
    V[ 4] = Point{3,Float64}(-0.5,  0.0,  -b)
    V[ 5] = Point{3,Float64}(-0.0, -0.5,   b)
    V[ 6] = Point{3,Float64}( 0.5, -0.0,   b)
    V[ 7] = Point{3,Float64}( 0.0,  0.5,   b)
    V[ 8] = Point{3,Float64}(-0.5,  0.0,   b)
    V[ 9] = Point{3,Float64}(-0.5, -0.5, 0.0)
    V[10] = Point{3,Float64}( 0.5, -0.5, 0.0)
    V[11] = Point{3,Float64}( 0.5,  0.5, 0.0)
    V[12] = Point{3,Float64}(-0.5,  0.5, 0.0)
    V[13] = Point{3,Float64}( 0.0,  0.0,  -a)
    V[14] = Point{3,Float64}( 0.0,  0.0,   a)

    F = Vector{QuadFace{Int}}(undef,12)
    F[ 1] = QuadFace{Int}( 1, 10, 5,  9)
    F[ 2] = QuadFace{Int}( 2, 11, 6, 10)
    F[ 3] = QuadFace{Int}( 3, 12, 7, 11)
    F[ 4] = QuadFace{Int}( 4,  9, 8, 12)
    F[ 5] = QuadFace{Int}( 5, 10, 6, 14)
    F[ 6] = QuadFace{Int}( 6, 11, 7, 14)
    F[ 7] = QuadFace{Int}( 7, 12, 8, 14)
    F[ 8] = QuadFace{Int}( 8,  9, 5, 14)
    F[ 9] = QuadFace{Int}( 1, 13, 2, 10)
    F[10] = QuadFace{Int}( 2, 13, 3, 11)
    F[11] = QuadFace{Int}( 3, 13, 4, 12)
    F[12] = QuadFace{Int}( 4, 13, 1,  9)

    if !isone(w)
        V .*= w
    end

    return F,V
end

"""
    tri2quad(F,V; method=:split)

Converts triangles to quads

# Description
This function converts the input triangular mesh, defined by the faces `F` and 
vertices `V`, to a quadrangulation. The method for this conversion is set using
the attribute `method` which can be set to `:split`, splitting each triangle 
into 3 quads by introducing a new central node, or `:rhombic`. whereby each 
triangle edge is used to construct a rhombic quadrilateral face. 
"""
function tri2quad(F,V; method=:split)
    # Get mesh edges 
    E = meshedges(F) # Non-unique edges
    Eu,_,inv_E = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
    Vf = simplexcenter(F,V) # Mid face points

    if method == :split
        # Compute new vertices 
        Ve = simplexcenter(Eu,V) # Mid edge points
        Vq = [V; Vf; Ve] # Append point sets

        # Create quadrilateral faces 
        Fq = Vector{QuadFace{Int}}(undef,length(F)*3)
        inv_E .+= length(V)+length(F) # Offset 
        nf = length(F)
        for (i,f) in enumerate(F)
            i_f = 1 + (i-1)*3 # Index of first face
            i_vf = i + length(V) # Index of mid face points
            Fq[i_f  ] = QuadFace{Int}(f[1], inv_E[i     ], i_vf, inv_E[i+2*nf] )
            Fq[i_f+1] = QuadFace{Int}(f[2], inv_E[i+  nf], i_vf, inv_E[i     ] )
            Fq[i_f+2] = QuadFace{Int}(f[3], inv_E[i+2*nf], i_vf, inv_E[i+  nf] )
        end
    elseif method == :rhombic
        con_E2F = con_edge_face(F,Eu,inv_E)
        nv = length(V)
        Fq = Vector{QuadFace{Int}}()
        seen = Dict{Int,LineFace{Int}}()  
        for (i,e) in enumerate(Eu)  
            ind_F = con_E2F[i]
            f1 = F[ind_F[1]]
            if length(ind_F) == 2 # Edge touches two faces 
                if e == f1[1:2] || e == f1[2:3] || e == f1[[3,1]] 
                    push!(Fq, QuadFace{Int}(e[1],ind_F[2]+nv,e[2],ind_F[1]+nv) )
                else
                    push!(Fq, QuadFace{Int}(e[1],ind_F[1]+nv,e[2],ind_F[2]+nv) )
                end 
            else # Edge touches one face
                i_f = ind_F[1] 
                if haskey(seen,i_f) # Seen before
                    ep = seen[i_f]
                    B = [in(j,e) && in(j,ep) for j in f1]
                    Vf[i_f] = V[f1[B][1]]
                else
                    Vf[i_f] = 0.5*(V[e[1]]+V[e[2]])
                    seen[i_f] = e
                end
            end    
        end
        Vq = [V; Vf] # Append point sets
    end
    return Fq,Vq
end

"""
    tetgenmesh(F::Array{NgonFace{N,TF}, 1},V::Vector{Point{3,TV}}; facetmarkerlist=nothing, V_regions=nothing,region_vol=nothing,V_holes=nothing,stringOpt="paAqYQ")  where N where TF<:Integer where TV<:Real

Creates a tetrahedral mesh

# Description
This function uses the TetGen.jl library to mesh the input geometry defined by 
the faces `F` and the vertices `V` using tetrahedral elements. Several optional
input parameters are available: 

* `facetmarkerlist`, a vector of integers with the same length as `F` and defines a face label for each face. 
* `V_regions`, a vector of points inside regions which require tetrahedral meshing.
* `region_vol`, a vector of scalar values to denote the desired tetrahedral volume for each region.  
* `V_holes`, a vector of points inside holes (voids) that should remain empty. 
* `stringOpt`, the TetGen command string to use. See also the TetGen documentation. 
"""
function tetgenmesh(F::Array{NgonFace{N,TF}, 1},V::Vector{Point{3,TV}}; facetmarkerlist=nothing, V_regions=nothing,region_vol=nothing,V_holes=nothing,stringOpt="paAqYQ")  where N where TF<:Integer where TV<:Real

    # Initialise TetGen input
    input = TetGen.RawTetGenIO{Cdouble}()

    # Add points     
    input.pointlist = reduce(hcat,V)

    # Add faces (facets)
    TetGen.facetlist!(input,F)

    # Add face boundary markers 
    if isnothing(facetmarkerlist)
        input.facetmarkerlist = ones(length(F))
    elseif length(facetmarkerlist) == length(F)
        input.facetmarkerlist = facetmarkerlist
    else
        throw(ArgumentError("The length of the facetmarkerlist should be length(F)"))        
    end

    # Add region definitions
    if isnothing(V_regions)
        V_regions = [mean(V)]
    end

    if isnothing(region_vol)        
        vol = mean(edgelengths(F,V))^3 / (6.0*sqrt(2.0)) # Volume for regular tetrahedron with mean edge length
        region_vol = fill(vol,length(V_regions))        
    elseif length(region_vol) == 1
        region_vol = fill(region_vol,length(V_regions))
    elseif length(region_vol) != length(V_regions)    
        throw(ArgumentError("The length of region_vol should be 1 or length(V_regions)"))
    end

    input.regionlist = Matrix{Cdouble}(undef,5,length(V_regions))
    for (i,v) in enumerate(V_regions)
        input.regionlist[:,i] = [v[1]; v[2]; v[3]; i; region_vol[i]]
    end
    
    # Add hole specifications 
    if !isnothing(V_holes)
        input.holelist = reduce(hcat,V_holes)
    end

    # TetGen meshing
    TET = tetrahedralize(input, stringOpt)

    # Retrieve elements, nodes, faces, and element/face labels
    E = [Tet4{TF}(e) for e in eachcol(TET.tetrahedronlist)] # Tetrahedral elements
    CE = vec(TET.tetrahedronattributelist)
    V = [Point{3,TV}(v) for v in eachcol(TET.pointlist)] # Vertices
    Fb = [TriangleFace{TF}(reverse(f)) for f in eachcol(TET.trifacelist)]
    Cb = TET.trifacemarkerlist

    return E,V,CE,Fb,Cb
end

"""
    surfacevolume(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real

Computes closed surface volume

# Description
This function computes the volume of a closed surface defined by the faces `F` 
and the vertices `V`. 
"""
function surfacevolume(F::Vector{NgonFace{N,TF}},V::Vector{Point{ND,TV}}) where N where TF<:Integer where ND where TV<:Real
    Z = [v[3] for v in V] # Z-coordinate for all vertices
    vol = 0.0 
    for f in F 
        zm = mean(Z[f])       
        @inbounds for qe in 1:N # Loop from first to end-1            
            # vol += zm * cross(V[f[qe]],V[f[mod1(qe+1,N)]])[3] # Add next edge contribution
            vol += zm * ( V[f[qe]][1]*V[f[mod1(qe+1,N)]][2] - V[f[qe]][2]*V[f[mod1(qe+1,N)]][1]) # Add next edge contribution          
        end 
    end    
    return vol/2.0
end

"""
    tetvolume(E::Vector{Tet4{T}},V::Vector{Point{ND,TV}}) where T<:Integer where ND where TV<:Real

Computes tetrahedral volumes

# Description 
This function computes the volume for each tetrahedron defined by the input `E`, 
a vector of Tet4 elements, and `V` the point coordinates. 
"""
function tetvolume(E::Vector{Tet4{T}},V::Vector{Point{ND,TV}}) where T<:Integer where ND where TV<:Real
    vol = Vector{TV}(undef,length(E))
    for (i,e) in enumerate(E)
        vol[i] = sum( (V[e[1]]-V[e[4]]) .* cross((V[e[3]]-V[e[4]]),(V[e[2]]-V[e[4]])) )
    end
    return vol/6.0
end

"""
    extrudefaces(F::Vector{NgonFace{NF,TF}},V::Vector{Point{ND,TV}}; extent=1.0, direction=:positive, num_steps=2, N::Union{Vector{Point{ND,TN}},Vector{Vec{ND, TN}},Nothing}=nothing) where NF where TF<:Integer where ND where TV<:Real where TN<:Real

Extrudes/thickens faces to form elements

# Description
The inputs surface mesh, defined by the faces `F` and vertices `V` is extruded to
create volumetric elements. Quadrilateral and triangular input faces are supported. 
These extrude into hexahedral and pentahedral elements respectively. 
The following input parameters are defined: 
* `extent<:Real` (default = 1.0) the length of the extrusion   
* `direction` is a symbol that is either `:positive` (default), `:negative`, or `:both`. 
* `N` The extrusion vectors. The default is nothing in which case the local 
vertex normals are used. 
* `num_steps` (default is 2) is the number of nodes in the extrude direction, the 
number of elements in the extrude direction is therefore `num_steps-1`. 
"""
function extrudefaces(F::Union{NgonFace{NF,TF},Vector{NgonFace{NF,TF}}},V::Vector{Point{ND,TV}}; extent=1.0, direction=:positive, num_steps=2, N::Union{Vector{Point{ND,TN}},Vector{Vec{ND, TN}},Nothing}=nothing) where NF where TF<:Integer where ND where TV<:Real where TN<:Real
    
    if isa(F,NgonFace{NF,TF})
        F = [F] 
    end

    # Compute normal directions if not provided
    if isnothing(N)
        N = vertexnormal(F,V; weighting=:area)
    end

    # Check if num_steps is okay
    if num_steps<1
        throw(ArgumentError("num_steps=$num_steps is not valid. num_steps should be larger than 0."))
    end

    # Set element type
    face_type=eltype(F)
    if face_type == QuadFace{TF}
        element_type = Hex8{TF}
    elseif face_type == TriangleFace{TF}
        element_type = Penta6{TF}
    else
        throw(ArgumentError("$face_type face type not supported. Supported types are QuadFace and TriangleFace."))
    end

    if direction == :positive # Default
        # no action needed
    elseif direction == :negative # Against N from V
        V -= extent .* N #Shift V in negative direction by half the extent      
    elseif direction == :both # Extrude both ways from V                
        V -= extent/2.0 .* N #Shift V in negative direction by half the extent      
    else
        throw(ArgumentError("$direction is not a valid direction, Use :positive, :negative, or :both.")) 
    end

    # Create coordinates and elements
    nv = length(V)
    nf = length(F)
    E = Vector{element_type}(undef,length(F)*(num_steps-1))
    VE = repeat(V,num_steps) # Vector{eltype(V)}(undef,n*m)
    for q = 1:num_steps-1
        # Create offset coordinates
        iv = 1 + q*nv
        VE[iv:(iv-1)+nv] += q/(num_steps-1)*N*extent

        # Create elements 
        ie = 1 + (q-1)*nf
        for (i,f) in enumerate(F)            
            E[ie+(i-1)] = (element_type)([f .+ (nv*(q-1)); f .+ (nv*q)])            
        end
    end

    return E, VE
end

"""
    filletcurve(V::Vector{Point{NV,TV}}; rMax::Union{Vector{T},T,Nothing}=nothing, constrain_method = :max, n=25, close_loop = false, eps_level = 1e-6) where TV<:Real where NV where T<:Real

Fillets/rounds curves

# Description
The function takes in a curve defined by the points `V` and applies filleting
(or rounding) to each "corner" (i.e. a point between two neighbouring points). 
The maximum radius `rMax` is used as the largest possible radius to use. If this, 
radius is not possible (e.g. if input points are too close), then a lower Radius
is used. 
"""
function filletcurve(V::Vector{Point{NV,TV}}; rMax::Union{Vector{T},T,Nothing}=nothing, constrain_method = :max, n=25, close_loop = false, eps_level = 1e-6) where TV<:Real where NV where T<:Real
    VP = deepcopy(V)
    if length(VP)>2
        if close_loop        
            VP = pushfirst!(VP,VP[end])
            VP = push!(VP,VP[2])   
            if isa(rMax,Vector{})
                rMax = pushfirst!(rMax,rMax[end])
                rMax = push!(rMax,rMax[2])   
            end
        end

        VC = Vector{eltype(VP)}()
        if !close_loop
            push!(VC,VP[1])
        end
        
        e1 = VP[2]-VP[1]
        e2 = VP[3]-VP[2]
        lp = 0.0
        i_last = length(VP)-1
        L = [norm(VP[i+1]-VP[i]) for i in 1:i_last]
        for i = 2:i_last
            if isa(rMax,Vector{})
                r = rMax[i]
            else 
                r = rMax
            end

            if i == 2
                e1 = VP[i]-VP[i-1]
            else
                e1 = e2
            end
            e2 = VP[i+1]-VP[i]        

            n1 = normalizevector(e1)
            n2 = normalizevector(e2)
            n3 = normalizevector(cross(n1,n2))
            α = pi - acos(dot(n1,n2))    
            if abs(α-pi) > eps_level        
                β = pi/2 - α/2.0

                n1_e = normalizevector(cross(n3,n1))
                n2_e = normalizevector(cross(n3,n2))
                m1 = normalizevector(n1_e + n2_e)
            
                l1 = L[i-1]
                l2 = L[i]     
                if constrain_method == :max
                    if i==2
                        if close_loop
                            lp = L[end]/2
                        else
                            lp=0.0
                        end
                    end
                    if i == i_last
                        if close_loop
                            l3 = L[2]
                        else
                            l3 = 0 
                        end
                    else
                        l3 = L[i+1] 
                    end
                                
                    if l3<l2   
                        l = min(l1-lp,l2-l3/2)
                    else
                        l = min(l1-lp,l2/2)    
                    end       
                elseif constrain_method == :mid
                    l = min(l1,l2)/2
                else
                    throw(ArgumentError("Incorrect constrain_method set $constrain_method, use :mid, or :max"))
                end

                rFit = l/tan(β)
                if !isnothing(r) && r<=rFit
                    rNow = r
                    lNow = rNow*tan(β) # Update as as radius used is smaller                                        
                    if isapprox(lNow,l,atol=eps_level)                        
                        fullRound = true
                    else
                        fullRound = false
                    end
                else
                    fullRound = true
                    rNow = rFit
                    lNow = l
                end

                if rNow > 0                    
                    d = abs(rNow/cos(β))
                    Q = RotMatrix3{Float64}(reduce(hcat,[-m1,normalizevector(cross(n3,-m1)),n3]))        
                                        
                    m1p = Q'*(-m1)
                    a = atan(m1p[2]/m1p[1])
                    
                    vc = (VP[i]+m1*d)                
                    Vc = [Q*Point{3, Float64}(rNow*cos(t),rNow*sin(t),0.0) + vc for t in range(a-β,a+β,n)]

                    if !close_loop && i==2 && fullRound && isapprox(norm(Vc[1]-VC[1]),0.0,atol=eps_level)                         
                        # Start is the same as first point of input curve
                        Vc = Vc[2:end] # Remove first point 
                    end
                    
                    if i>2 && fullRound && isapprox(norm(Vc[1]-VC[end]),0.0,atol=eps_level)                        
                        # The start of the current segment is the same as last point built curve so far 
                        Vc = Vc[2:end] # Remove first point 
                    end   

                    if close_loop && i == i_last && fullRound && isapprox(norm(Vc[end]-VC[1]),0.0,atol=eps_level)                        
                        # The end of the current segment is the same as the first point on the closed loop curve built
                        Vc = Vc[1:end-1] # Remove last point
                    end

                    append!(VC,Vc) # Append new curve segment to curve being constructed
                else
                    push!(VC,VP[i]) # No rounding so just add point
                    lNow=0.0                    
                end
                lp = lNow # Store as previous l 
            else
                push!(VC,VP[i]) # No rounding so add point
                lNow=0.0  
            end    
        end
        if !close_loop
            if isapprox(norm(VC[end]-VP[end]),0.0,atol=eps_level)                
                VC[end] = VP[end] # overwrite end
            else
                push!(VC,VP[end]) # add end
            end
        end
        return VC
    else
        return VP
    end
end

"""
    squircle(r::T,n::Int,τ=0.5; atol=1e-6, dir=:acw) where T <: Real

Creates the squircle curve

# Description
This function returns `n` points on a squircle. The squircle curve is defined 
using the radius `r`, and the parameter `τ`. The latter controls the morphing 
between a circle (`τ=0`) and square (`τ=1`).
"""
function squircle(r::T,n::Int,τ=0.5; atol=1e-6, dir=:acw) where T <: Real
    if isapprox(τ,0.0,atol=atol)
        V = circlepoints(r,n, dir=dir)
    else 
        p = circlerange(n;dir=dir)

        # Compute s from τ
        st = r^2*((sqrt(2.0)/2.0*(1.0-τ))+τ)^2
        s = (r/st)*sqrt(2.0*st-r^2)

        V = Vector{Point{3,Float64}}(undef,n)
        V[1] = [r,0.0,0.0]    
        for i in 2:1:n
            if isapprox(p[i],0.5*π,atol=atol)
                V[i] = [0.0,r,0.0]
            elseif isapprox(p[i],π,atol=atol)
                V[i] = [-r,0.0,0.0]
            elseif isapprox(p[i],1.5*π,atol=atol)
                V[i] = [0.0,-r,0.0]
            else
                cos_p = cos(p[i])
                sin_p = sin(p[i])                
                q = (r/(s*sqrt(2.0))) *  sqrt(1.0 - sqrt.(1.0 - s^2 * max(0.0,min(1.0,(2.0 * sin_p * cos_p)^2))))
                x = sign(cos_p) / abs(sin_p)*q
                y = sign(sin_p) / abs(cos_p)*q            
                V[i] = [x,y,0]
            end
        end        
    end
    return V
end

"""
    circlerange(n::Int; dir=:acw, deg=false)

Creates circular angles

# Description
This function returns `n` angles for an even circular distribution of points. 
The optional input `dir` can be set to `:acw` (default) resulting in an 
anti-clockwise set of angles, and can be set to `:cw` for a clockwise set of 
angles. Angles are returned in radians since `deg` is `false` by default. 
Using `deg=true` results in angles in degrees. 
"""
function circlerange(n::Int; dir=:acw, deg=false)
    if dir==:acw
        if deg 
            return range(0.0, 360 - 360/n, n)
        else
            return range(0.0, 2.0*π - (2.0*π)/n, n)
        end
    elseif dir==:cw
        if deg 
            return range(0.0, 360/n - 360, n)
        else
            return range(0.0, (2.0*π)/n - 2.0*π, n)
        end
    else
        throw(ArgumentError("Invalid dir specified :$dir, use :acw or :cw"))
    end
end


"""
    edgefaceangles(F::Vector{NgonFace{NF,TF}},V::Vector{Point{ND,TV}}; deg=false) where NF where TF<:Integer where ND where TV<:Real 

Computed angles between faces

# Description
This function computes the angle between two faces for each unique edge in the 
mesh specified by `F` and the vertices `V`. If the input mesh consists of n 
unique edges then the output features n angles. For boundary edges, where no
pair of faces exists, the angle returned is NaN. The default behaviour results
in angles being computed in radians. However by specifying `deg=true` the user
can request degrees instead. Finally two additional outputs are created, namely
the unique edge vector `E_uni` as well as the edge-to-face connectivity vector
`,con_E2F`. 
"""
function edgefaceangles(F::Vector{NgonFace{NF,TF}},V::Vector{Point{ND,TV}}; deg=false) where NF where TF<:Integer where ND where TV<:Real 
    E = meshedges(F) # The non-unique mesh edges 
    E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true) # Unique mesh edges and inverse indices
    con_E2F = con_edge_face(F,E_uni,indReverse) # The edge-to-face connectivity

    if length(F)>1 # More than one face so compute connectivity
        A = fill(TV(NaN),length(E_uni)) # Vector for storing angles for each edge
        for i_e in eachindex(E_uni) # Looping over all unique edges (e.g. 1-2 is the same as 2-1)
            if length(con_E2F[i_e])==2 # If the edge is not a boundary edge it touches 2 faces and angles can be computed           
                e = E_uni[i_e] # The current edge
                i_f1 = con_E2F[i_e][1] # The index of the first face sharing the current edge
                i1 = findfirst(F[i_f1].==e[1]) # The index of the first edge point
                i2 = findfirst(F[i_f1].==e[2]) # The index of the second edge point                           
                n1 = facenormal(F[con_E2F[i_e][1]],V) # The normal direction for the first face 
                n2 = facenormal(F[con_E2F[i_e][2]],V)  # The normal direction for the second face         
                ne = V[E_uni[i_e][2]]-V[E_uni[i_e][1]] # Current edge vector                     
                if i2==mod1(i1+1,length(F[1])) # If the second point is "next" for the current edge
                    s = sign(dot(cross(n2,n1),ne))    
                else # Inverse situation
                    s = sign(dot(cross(n1,n2),ne))
                end
                
                # Compute the face angles 
                if deg # Compute angles in degrees                       
                    A[i_e] = s*acosd(clamp(dot(n1,n2),-1.0,1.0))
                else # Compute angles in radians
                    A[i_e] = s*acos(clamp(dot(n1,n2),-1.0,1.0))                         
                end                 
            end
        end
        return A,E_uni,con_E2F
    else # Just one face, so return NaNs for the angles
        return fill(NaN,length(E)), E_uni,con_E2F
    end
end


"""
    faceanglesegment(F::Vector{NgonFace{NF,TF}},V::Vector{Point{ND,TV; deg=false, angleThreshold = pi/8, indStart = 1)  where NF where TF<:Integer where ND where TV<:Real 

Segments surfaces using face angles

# Description
This function takes in a surface mesh defined by the faces `F` and the vertices 
`V`, and segments the surface mesh based on face angles. The output consists of 
a "feature label" vector `G` (a Vector{Int}, with the same length of `F`) 
whereby adjacent faces whosee angle is smaller or equal to the `angleThreshold` 
(default is pi/8) receive the same label. Hence this function allows one to find 
all faces with a similar orientation, for instance all top or side faces of some 
mesh. The function uses radians by default. However, buy specifying the optional 
parameter `deg=true` the user request that angles and the `angleThreshold` are
in degrees.  
"""
function faceanglesegment(F::Vector{NgonFace{NF,TF}},V::Vector{Point{ND,TV}}; deg=false, angleThreshold = pi/8, indStart = 1)  where NF where TF<:Integer where ND where TV<:Real 

    A,_,con_E2F = edgefaceangles(F,V; deg=deg)
    indToCheck = [indStart]
    indToCheckNew = Vector{Int}()
    G = zeros(Int,length(F)) # Group labels
    G[indStart] = 1 
    g = 1

    while true 
        while true 
            for indFaceNow in indToCheck # Loop over all indices of faces to check on 
                while true # While we are finding unvisited edges attached to the current face
                    indEdgeFound = findfirst(x -> in(indFaceNow,x),con_E2F) # index of edge attached to face. check empty
                    if !isnothing(indEdgeFound) # If we found an edge 
                        if abs(A[indEdgeFound])<=angleThreshold # Is angle at this edge low enough
                            for indFace in con_E2F[indEdgeFound] # then add the other face
                                if indFace!=indFaceNow # If other face                                    
                                    push!(indToCheckNew,indFace) # Add to check for next round
                                    G[indFace] = g # Assign group label
                                end
                            end
                        end
                        con_E2F[indEdgeFound] = [] # Empty visited edge
                    else
                        break
                    end
                end
            end

            if isempty(indToCheckNew)
                g+=1
                indToCheck = findfirst(x -> x==0,G)
                break
            else
                indToCheck = deepcopy(indToCheckNew) # Create new set to check up on 
                indToCheckNew = Vector{Int}() # Create new set empty again
            end
        end        

        if isnothing(indToCheck)
            break
        else
            G[indToCheck] = g # Assign group label
        end
    end
    return G
end

"""
    eulerchar(F,V=nothing; E=nothing)

Computes the Euler characteristic

# Description
This function computes the Euler characteristic for the input surface defined by 
the faces `F` and vertices `V`. The edges `E` are on optional input. 
The Euler characteristic is defined as: 
`X = nV-nE-nF`
, where `nV`, `nE`, and `nF` define the number of surface vertices, edges, and 
faces respectively. It is assumed all inputs are set of unique entities, e.g. 
no vertices, edges, or faces occur multiple times. 
"""
function eulerchar(F,V=nothing,E=nothing)
    nf = length(F)
    if isnothing(V)
        nv = maximum(reduce(vcat,F)) # Use largest index (assumes all points used in mesh)
    else
        nv = length(V)
    end
    if isnothing(E)
        ne = length(meshedges(F; unique_only=true))        
    else 
        ne = length(E)
    end    
    return nv-ne+nf
end

"""
    rhombicdodecahedronfoam(w::T,n::Union{Tuple{Vararg{Int, 3}}, Array{Int, 3}}; merge = true, orientation = :allign) where T<:Real

Creates a rhombicdodecahedron foam

# Description
This function creates a rhombicdodecahedron foam structure with a desired number
of cells in each direction. The input is the cell width `w` and a 1-by-3 tuple 
`n` defining the number of cells in each direction. The output consists of the 
set of rhombic dodecahedron elements `E` and their vertex coordinates `V`. 
"""
function rhombicdodecahedronfoam(w::T,n::Union{Tuple{Vararg{Int, 3}}, Array{Int, 3}}; merge = true, orientation = :allign) where T<:Real

    if orientation == :allign
        # Create vertices of single rhombic dodecahedron with "radius" 1
        a = sqrt(2)/2
        b = sqrt(2)/4
        
        V0 = Vector{Point{3,Float64}}(undef,14)
        V0[ 1] = Point{3,Float64}(-0.0, -0.5,  -b)
        V0[ 2] = Point{3,Float64}( 0.5, -0.0,  -b)
        V0[ 3] = Point{3,Float64}( 0.0,  0.5,  -b)
        V0[ 4] = Point{3,Float64}(-0.5,  0.0,  -b)
        V0[ 5] = Point{3,Float64}(-0.0, -0.5,   b)
        V0[ 6] = Point{3,Float64}( 0.5, -0.0,   b)
        V0[ 7] = Point{3,Float64}( 0.0,  0.5,   b)
        V0[ 8] = Point{3,Float64}(-0.5,  0.0,   b)
        V0[ 9] = Point{3,Float64}(-0.5, -0.5, 0.0)
        V0[10] = Point{3,Float64}( 0.5, -0.5, 0.0)
        V0[11] = Point{3,Float64}( 0.5,  0.5, 0.0)
        V0[12] = Point{3,Float64}(-0.5,  0.5, 0.0)
        V0[13] = Point{3,Float64}( 0.0,  0.0,  -a)
        V0[14] = Point{3,Float64}( 0.0,  0.0,   a)

        nv = length(V0) # Number of vertices

        # Scale to obtain desired radius if needed
        if !isone(w)
            V0.*=w
        end

        # Create rhombicdodecahedron element set    
        e0 = Rhombicdodeca14{Int}(1:nv)
        m = ceil(Int,n[3]/2)    
        k = m*(n[1]*n[2]) + (n[3]-m)*((n[1]-1)*(n[2]-1))
        E = Vector{Rhombicdodeca14{Int}}(undef,k)    
        @inbounds for i in 1:k        
            E[i] = e0.+(i-1)*nv
        end       

        # Create coordinates 
        xRange = range(0.0,(n[1]-1)*w,n[1])
        yRange = range(0.0,(n[2]-1)*w,n[2])
        zRange = range(0.0,(n[3]-1)*sqrt(2)*w/2,n[3])    
        V = Vector{Point{3,Float64}}(undef,k*nv)

        i = 1
        for (i_z,z) in enumerate(zRange)
            if iseven(i_z)                  
                xRange_now = xRange[1:end-1] .+ w/2
                yRange_now = yRange[1:end-1] .+ w/2            
            else
                xRange_now = xRange
                yRange_now = yRange
            end

            for x = xRange_now
                for y = yRange_now                          
                    V[i:(i+nv-1)] = [v+Point{3,Float64}(x,y,z) for v in V0]                
                    i += nv
                end
            end        
        end
    elseif orientation == :diagonal
        # Create vertices of single rhombic dodecahedron with width 
        V0 = Vector{Point{3,Float64}}(undef,14)
        V0[ 1] = Point{3,Float64}(-0.25, -0.25, -0.25) 
        V0[ 2] = Point{3,Float64}( 0.25, -0.25, -0.25) 
        V0[ 3] = Point{3,Float64}( 0.25,  0.25, -0.25) 
        V0[ 4] = Point{3,Float64}(-0.25,  0.25, -0.25) 
        V0[ 5] = Point{3,Float64}(-0.25, -0.25,  0.25) 
        V0[ 6] = Point{3,Float64}( 0.25, -0.25,  0.25) 
        V0[ 7] = Point{3,Float64}( 0.25,  0.25,  0.25) 
        V0[ 8] = Point{3,Float64}(-0.25,  0.25,  0.25) 
        V0[ 9] = Point{3,Float64}(-0.50, -0.00,  0.00)
        V0[10] = Point{3,Float64}( 0.00, -0.50,  0.00)
        V0[11] = Point{3,Float64}( 0.50,  0.00,  0.00)
        V0[12] = Point{3,Float64}( 0.00,  0.50,  0.00)
        V0[13] = Point{3,Float64}( 0.00, -0.00, -0.50)
        V0[14] = Point{3,Float64}( 0.00, -0.00,  0.50)

        # Scale to obtain desired radius if needed
        if !isone(w)
            V0.*=w
        end

        # Create rhombicdodecahedron element set    
        e0 = Rhombicdodeca14{Int}(1:14)
        k = prod(n) + (n[1]-1) * (n[2]-1) * n[3] + (n[1]-1) * n[2] * (n[3]-1) + n[1] * (n[2]-1)  * (n[3]-1)
        E = Vector{Rhombicdodeca14{Int}}(undef,k)    
        @inbounds for i in 1:k        
            E[i] = e0.+(i-1)*14
        end       

        # Create coordinates 
        x = range(0.0,(n[1]-1)*w,n[1])
        y = range(0.0,(n[2]-1)*w,n[2])
        z = range(0.0,(n[3]-1)*w,n[3])    
        V = Vector{Point{3,Float64}}(undef,k*14)
        i = 1
        for x = x
            for y = y
                for z = z
                    V[i:(i+14-1)] = [v+Point{3,Float64}(x,y,z) for v in V0]                
                    i += 14
                end
            end
        end
        for x = x[1:end-1].+w/2
            for y = y
                for z = z[1:end-1].+w/2
                    V[i:(i+14-1)] = [v+Point{3,Float64}(x,y,z) for v in V0]                
                    i += 14
                end
            end
        end
        for x = x
            for y = y[1:end-1].+w/2
                for z = z[1:end-1].+w/2
                    V[i:(i+14-1)] = [v+Point{3,Float64}(x,y,z) for v in V0]                
                    i += 14
                end
            end
        end
        for x = x[1:end-1].+w/2
            for y = y[1:end-1].+w/2
                for z = z
                    V[i:(i+14-1)] = [v+Point{3,Float64}(x,y,z) for v in V0]                
                    i += 14
                end
            end
        end
    end
    
    # merge vertices if requested
    if merge
        V,_,indMap = mergevertices(V; pointSpacing=w)
        indexmap!(E,indMap)        
    end
    return E,V
end

"""
    truncatedoctahedron(w=1.0)

Creates a truncated octahedron

# Description
This function creates a truncated octahedron. The input cell width `w` is used 
to define the cell faces `F` and vertices `V`. 
"""
function truncatedoctahedron(w=1.0)    
    # Coordinates
    V = Vector{Point{3,Float64}}(undef,24)
    V[ 1] = Point{3,Float64}( 0.50, -0.25,  0.00) 
    V[ 2] = Point{3,Float64}(-0.25, -0.50,  0.00) 
    V[ 3] = Point{3,Float64}(-0.00, -0.50, -0.25) 
    V[ 4] = Point{3,Float64}(-0.00, -0.25,  0.50) 
    V[ 5] = Point{3,Float64}( 0.25,  0.50,  0.00) 
    V[ 6] = Point{3,Float64}( 0.25, -0.00, -0.50) 
    V[ 7] = Point{3,Float64}( 0.50, -0.00,  0.25) 
    V[ 8] = Point{3,Float64}(-0.50,  0.25,  0.00) 
    V[ 9] = Point{3,Float64}( 0.00,  0.25, -0.50) 
    V[10] = Point{3,Float64}( 0.00,  0.50,  0.25) 
    V[11] = Point{3,Float64}(-0.25,  0.00, -0.50) 
    V[12] = Point{3,Float64}(-0.50,  0.00,  0.25) 
    V[13] = Point{3,Float64}( 0.25, -0.50,  0.00) 
    V[14] = Point{3,Float64}(-0.50, -0.25,  0.00) 
    V[15] = Point{3,Float64}(-0.00, -0.25, -0.50) 
    V[16] = Point{3,Float64}(-0.00, -0.50,  0.25) 
    V[17] = Point{3,Float64}( 0.50,  0.25,  0.00) 
    V[18] = Point{3,Float64}( 0.50, -0.00, -0.25) 
    V[19] = Point{3,Float64}( 0.25, -0.00,  0.50) 
    V[20] = Point{3,Float64}(-0.25,  0.50,  0.00) 
    V[21] = Point{3,Float64}( 0.00,  0.50, -0.25) 
    V[22] = Point{3,Float64}( 0.00,  0.25,  0.50) 
    V[23] = Point{3,Float64}(-0.50,  0.00, -0.25) 
    V[24] = Point{3,Float64}(-0.25,  0.00,  0.50)

    # Scale width if needed
    if !isone(w)
        V .*= w
    end

    # Hexagonal faces
    F1 = Vector{NgonFace{6,Int}}(undef,8)
    F1[ 1] = NgonFace{6,Int}( 6, 18,  1, 13,  3, 15)
    F1[ 2] = NgonFace{6,Int}( 9, 21,  5, 17, 18,  6)
    F1[ 3] = NgonFace{6,Int}(11, 23,  8, 20, 21,  9)
    F1[ 4] = NgonFace{6,Int}(15,  3,  2, 14, 23, 11)
    F1[ 5] = NgonFace{6,Int}(13,  1,  7, 19,  4, 16)
    F1[ 6] = NgonFace{6,Int}(17,  5, 10, 22, 19,  7)
    F1[ 7] = NgonFace{6,Int}(20,  8, 12, 24, 22, 10)
    F1[ 8] = NgonFace{6,Int}(14,  2, 16,  4, 24, 12)

    # Quadrilateral faces
    F2 = Vector{QuadFace{Int}}(undef,6)
    F2[ 1] = QuadFace{Int}(2,  3, 13, 16)
    F2[ 2] = QuadFace{Int}(17,  7,  1, 18)
    F2[ 3] = QuadFace{Int}(20, 10,  5, 21)
    F2[ 4] = QuadFace{Int}(14, 12,  8, 23)
    F2[ 5] = QuadFace{Int}(11,  9,  6, 15)
    F2[ 6] = QuadFace{Int}(4, 19, 22, 24)

    # Collect faces in tuple
    F = (F1,F2)

    return F,V
end

"""
    kelvinfoam(w::T,n::Union{Tuple{Vararg{Int, 3}}, Array{Int, 3}}; merge = true) where T<:Real

Creates a Kelvin foam

# Description
This function creates a Kelvin foam structure with a desired number of cells in
each direction. The input is the cell width `w` and a 1-by-3 tuple `n` defining
the number of cells in each direction. The output consists of the set of 
truncated octahedron elements `E` and their vertex coordinates `V`. 
"""
function kelvinfoam(w::T,n::Union{Tuple{Vararg{Int, 3}}, Array{Int, 3}}; merge = true) where T<:Real

    # Create vertices of single rhombic dodecahedron with "radius" 1    
    V0 = Vector{Point{3,Float64}}(undef,24)
    V0[ 1] = Point{3,Float64}( 0.50, -0.25,  0.00) 
    V0[ 2] = Point{3,Float64}(-0.25, -0.50,  0.00) 
    V0[ 3] = Point{3,Float64}(-0.00, -0.50, -0.25) 
    V0[ 4] = Point{3,Float64}(-0.00, -0.25,  0.50) 
    V0[ 5] = Point{3,Float64}( 0.25,  0.50,  0.00) 
    V0[ 6] = Point{3,Float64}( 0.25, -0.00, -0.50) 
    V0[ 7] = Point{3,Float64}( 0.50, -0.00,  0.25) 
    V0[ 8] = Point{3,Float64}(-0.50,  0.25,  0.00) 
    V0[ 9] = Point{3,Float64}( 0.00,  0.25, -0.50) 
    V0[10] = Point{3,Float64}( 0.00,  0.50,  0.25) 
    V0[11] = Point{3,Float64}(-0.25,  0.00, -0.50) 
    V0[12] = Point{3,Float64}(-0.50,  0.00,  0.25) 
    V0[13] = Point{3,Float64}( 0.25, -0.50,  0.00) 
    V0[14] = Point{3,Float64}(-0.50, -0.25,  0.00) 
    V0[15] = Point{3,Float64}(-0.00, -0.25, -0.50) 
    V0[16] = Point{3,Float64}(-0.00, -0.50,  0.25) 
    V0[17] = Point{3,Float64}( 0.50,  0.25,  0.00) 
    V0[18] = Point{3,Float64}( 0.50, -0.00, -0.25) 
    V0[19] = Point{3,Float64}( 0.25, -0.00,  0.50) 
    V0[20] = Point{3,Float64}(-0.25,  0.50,  0.00) 
    V0[21] = Point{3,Float64}( 0.00,  0.50, -0.25) 
    V0[22] = Point{3,Float64}( 0.00,  0.25,  0.50) 
    V0[23] = Point{3,Float64}(-0.50,  0.00, -0.25) 
    V0[24] = Point{3,Float64}(-0.25,  0.00,  0.50)

    nv = length(V0) # Number of vertices

    # Scale to obtain desired radius if needed
    if !isone(w)
        V0.*=w
    end

    # Create truncated octahedron element set    
    e0 = Truncatedocta24{Int}(1:nv)
    m = ceil(Int,n[3]/2)    
    k = m*(n[1]*n[2]) + (n[3]-m)*((n[1]-1)*(n[2]-1))
    E = Vector{Truncatedocta24{Int}}(undef,k)    
    @inbounds for i in 1:k        
        E[i] = e0.+(i-1)*nv
    end       

    # Create coordinates 
    xRange = range(0.0,(n[1]-1)*w,n[1])
    yRange = range(0.0,(n[2]-1)*w,n[2])
    zRange = range(0.0,(n[3]-1)*w/2,n[3])    
    V = Vector{Point{3,Float64}}(undef,k*nv)

    i = 1
    for (i_z,z) in enumerate(zRange)
        if iseven(i_z)                  
            xRange_now = xRange[1:end-1] .+ w/2
            yRange_now = yRange[1:end-1] .+ w/2            
        else
            xRange_now = xRange
            yRange_now = yRange
        end

        for x = xRange_now
            for y = yRange_now                          
                V[i:(i+nv-1)] = [v+Point{3,Float64}(x,y,z) for v in V0]                
                i += nv
            end
        end        
    end
      
    # merge vertices if requested
    if merge 
        V,_,indMap = mergevertices(V; pointSpacing=w)
        indexmap!(E,indMap)        
    end
    return E,V
end

"""
    minp(V::Vector{Point{N,T}}) where N where T <:Real   

Returns minimum coordinates

# Description
This function computes the minimum coordinates for all points. Points can be 
N-dimensional and the output is another point of the same dimensionality but 
with the lowest coordinate value for each direction. 
"""
function minp(V::Vector{Point{N,T}}) where N where T <:Real   
    m = fill(T(Inf),N)
    @inbounds for v in V
        @inbounds for j = 1:N            
            m[j] = min(m[j],v[j])
        end        
    end    
    return Point{N,T}(m)
end

"""
    maxp(V::Vector{Point{N,T}}) where N where T <:Real   

Returns maximum coordinates

# Description
This function computes the maximum coordinates for all points. Points can be 
N-dimensional and the output is another point of the same dimensionality but 
with the highest coordinate value for each direction. 
"""
function maxp(V::Vector{Point{N,T}}) where N where T <:Real   
    m = fill(T(-Inf),N)
    @inbounds for v in V
        @inbounds for j = 1:N            
            m[j] = max(m[j],v[j])
        end        
    end    
    return Point{N,T}(m)
end

"""
    ntrapezohedron(n,r=1.0)
    
Constructs an n-trapezohedron

# Description
This function creates the faces `F` and vertices `V` for an n-trapezohedron. 

The implementation is based on the equations presented here: 
https://mathworld.wolfram.com/Trapezohedron.html
"""
function ntrapezohedron(n,r=1.0)    
    m = 2*n # Number of faces
    R = 0.5*csc(pi/n)*r # Radius scaled to obtain desired midradius
    zc = r * sqrt(4-sec(pi/m)^2)/(4+8*cos(pi/n)) # z-offset of equatorial points
    zt = r * 1/4*cos(pi/m)*cot(pi/m)*csc(3*pi/m)*sqrt(4-sec(pi/m)^2) # z-offset of poles

    # Create equatorial points
    T = circlerange(m; dir=:acw) # Range of angles for points in circle
    V = Vector{Point{3,Float64}}(undef,m+2)
    for (i,t) in enumerate(T)
        if iseven(i)
            s = 1.0
        else
            s = -1.0
        end
        V[i] = Point{3, Float64}(R*cos(t),R*sin(t),s*zc) 
    end

    # Create pole points
    V[end-1] = Point{3, Float64}(0.0,0.0,-zt) 
    V[end  ] = Point{3, Float64}(0.0,0.0, zt) 

    # Construct faces
    F =  Vector{QuadFace{Int}}(undef,m)
    for i = 1:n # Use n steps to create m=2*n faces
        j = 1 + 2*(i-1) # First point index for bottom quads       
        F[i  ] = QuadFace{Int}(mod1(j+2,m),         j+1,           j, m+1)  # Bottom quad 
        F[i+n] = QuadFace{Int}(        j+1, mod1(j+2,m), mod1(j+3,m), m+2) # Top quad   
    end
    
    return F, V
end

"""
    spacing2numvertices(F::Vector{TriangleFace{TF}},V::Vector{Point{ND,TV}},pointSpacing::TP) where TF<:Integer where ND where TV<:Real where TP<:Real 

Point numbers from spacing

# Description 
This function helps to determine what number of vertices to resample a surface 
with to obtain a desired point spacing. The input consists of in initial surface
, defined by the faces `F` and the vertices `V`, and also the desired point 
spacing `pointSpacing`. Next the function uses the surface Euler characteristic
as well as knowledge of face area and point spacing changes for homogeneous 
face splitting (e.g. via `subtri``), to determine the theoretical number of 
points `NV` a resampled surface should have to present with the desired point
spacing. 
"""
function spacing2numvertices(F::Vector{TriangleFace{TF}},V::Vector{Point{ND,TV}},pointSpacing::TP) where TF<:Integer where ND where TV<:Real where TP<:Real 
    # The following assumes subtri like splitting whereby each edge spawns new 
    # point and each face is split into 4 new faces.  

    # Compute desired number of faces using area
    A_total = sum(facearea(F,V)) # Total surface area estimation
    A_ideal = (pointSpacing.^2.0*sqrt(3.0))/4.0 # Theoretical area of equilateral triangle
    nf_desired = (A_total/A_ideal) # Desired number of faces
    
    # Determine theoretical number of refinement iterations to use to get the desired number of faces
    E = meshedges(F; unique_only = true) # The unique mesh edges
    nf = length(F) # Number of input faces
    ne = length(E) # Number of input edges
    nRefineScalar = (log(nf_desired)-log(nf))/log(4) 

    if nRefineScalar<0    
        n = 0:-1:floor(nRefineScalar)-1
    else    
        n = 0:1:ceil(nRefineScalar)+1
    end
    
    # Use knowledge of triangle splitting to get corresponding sets of edge and face counts
    nf_sim = nf .* 4.0 .^n    
    ne_sim = zeros(length(n))
    ne_sim[1] = ne
    if nRefineScalar>0.0
        for i in 2:1:length(n)
            ne_sim[i] = 2.0 .* ne_sim[i-1] .+ (3.0*nf_sim[i-1])
        end
    else
        for i in 2:1:length(n)
            ne_sim[i] = 0.5 .* ( ne_sim[i-1] .- (3.0*nf_sim[i]) )
        end
    end

    # Use Euler's characteristic to determine theoretical number of vertices
    nv_sim = eulerchar(F,V) .+ ne_sim .- nf_sim

    
    # Interpolate at desired number of faces to get desired number of vertices 
    if nRefineScalar>0.0
        nv_Out = ceil(Int,lerp(n,nv_sim,nRefineScalar))
    else
        nv_Out = ceil(Int,lerp(reverse(collect(n)),reverse(nv_sim),nRefineScalar))
    end
    return nv_Out 
end

"""
    joingeom(G...)

Joins geometry

# Description 
This function joins geometry defined for instance by multiple face and vertex 
sets into one such set. All geometry such be of the same type such that they can
be joined. The input can for instance be n-sets of faces (or elements) and 
vertices e.g. appearing as inpus as: `F1,V1,F2,V2,...,FN,VN`.  
"""
function joingeom(G...)
    # The input G is of the form F1,V1,F2,V2,...FN,VN    
    if iseven(length(G)) # G should be even in length                         
        if length(G)==2
            # Only one set of faces and vertices
            return G[1],G[2],ones(length(G[1]))
        else
            # Multiple sets of faces and vertices
            F = deepcopy(G[1]) # Initial face set 
            V = deepcopy(G[2]) # Initial vertex set 
            C = ones(length(G[1]))                            
            indexShift = length(V)         
            c = 1   
            @inbounds for i in 3:1:length(G)
                g = G[i]
                if iseven(i)
                    append!(V,g) 
                    indexShift += length(g)
                else
                    c += 1
                    append!(F,[f.+indexShift for f in g])                     
                    append!(C,fill(c,length(g)))                                         
                end
            end
            return F,V,C        
        end
    else        
        throw(ArgumentError("Number of element sets does not seem to match the number of vertex sets"))        
    end
end

"""
    quadbox(boxDim,boxEl)

Creates quadrilateral box mesh

# Description 
This function uses the dimensions defined in `boxDim`, and the number of 
elements listed for each direction in `boxEl` to create a quadrilateral mesh for 
a box. The output consists of the faces `F`, the vertices `V`, and a face 
labelling `C` (for the 6 sides of the box). 
"""
function quadbox(boxDim,boxEl)
    F12,V12 = quadplate(boxDim[[1,2]],boxEl[[1,2]]; orientation=:up)
    F22,V22 = quadplate(boxDim[[1,3]],boxEl[[1,3]]; orientation=:up)
    F32,V32 = quadplate(boxDim[[3,2]],boxEl[[3,2]]; orientation=:up)
    return _faces2box(F12,V12,F22,V22,F32,V32,boxDim)
end

"""
    tribox(boxDim,pointSpacing)

Creates triangulated box mesh

# Description 
This function uses the dimensions defined in `boxDim`, and the desired point
spacing defined by `pointSpacing`, to create a triangulated mesh for a box. 
The output consists of the faces `F`, the vertices `V`, and a face 
labelling `C` (for the 6 sides of the box). 
"""
function tribox(boxDim,pointSpacing)
    xr = range(-boxDim[1]/2.0,boxDim[1]/2.0,1+ceil(Int64,boxDim[1]/pointSpacing))
    yr = range(-boxDim[2]/2.0,boxDim[2]/2.0,1+ceil(Int64,boxDim[2]/pointSpacing))
    zr = range(-boxDim[3]/2.0,boxDim[3]/2.0,1+ceil(Int64,boxDim[3]/pointSpacing))

    Vc = [Point{3,Float64}(x,-boxDim[2]/2,0.0) for x in xr]
    v1 = [Point{3,Float64}(boxDim[1]/2,y,0.0) for y in yr]
    v2 = [Point{3,Float64}(x, boxDim[2]/2,0.0) for x in xr]
    v3 = [Point{3,Float64}(-boxDim[1]/2,y,0.0) for y in yr]
    append!(Vc,v1[2:end])
    append!(Vc,reverse(v2[1:end-1]))
    append!(Vc,reverse(v3[2:end-1]))
    F12,V12,_ = regiontrimesh((Vc,),([1],),(pointSpacing))
    
    Vc = [Point{3,Float64}(x,-boxDim[3]/2,0.0) for x in xr]
    v1 = [Point{3,Float64}(boxDim[1]/2,z,0.0) for z in zr]
    v2 = [Point{3,Float64}(x, boxDim[3]/2,0.0) for x in xr]
    v3 = [Point{3,Float64}(-boxDim[1]/2,z,0.0) for z in zr]
    append!(Vc,v1[2:end])
    append!(Vc,reverse(v2[1:end-1]))
    append!(Vc,reverse(v3[2:end-1]))
    F22,V22,_ = regiontrimesh((Vc,),([1],),(pointSpacing))

    Vc = [Point{3,Float64}(z,-boxDim[2]/2,0.0) for z in zr]
    v1 = [Point{3,Float64}(boxDim[3]/2,y,0.0) for y in yr]
    v2 = [Point{3,Float64}(z,boxDim[2]/2,0.0) for z in zr]
    v3 = [Point{3,Float64}(-boxDim[3]/2,y,0.0) for y in yr]
    append!(Vc,v1[2:end])
    append!(Vc,reverse(v2[1:end-1]))
    append!(Vc,reverse(v3[2:end-1]))
    F32,V32,_ = regiontrimesh((Vc,),([1],),(pointSpacing))

    return _faces2box(F12,V12,F22,V22,F32,V32,boxDim)
end

"""
    _faces2box(F12,V12,F22,V22,F32,V32,boxDim)

Converts face set to box

# Description 
This unexported function helps turn 3 sets of faces into a mesh for a box. It is 
assumed that the faces are for the sides of a box and that each set is defined 
in the XY plane centered around 0,0. Hence to form the box these faces are 
rotated and shifted (and inverted if needed). Next the vertices and faces are
joined and merged to create a "watertight" vertex-sharing closed box. 
"""
function _faces2box(F12,V12,F22,V22,F32,V32,boxDim)
    V11 = V12 .+ Point{3, Float64}([0.0,0.0,-boxDim[3]/2])
    V12 = V12 .+ Point{3, Float64}([0.0,0.0, boxDim[3]/2])
    F11 = invert_faces(F12)

    Q = RotXYZ(-0.5*pi,0.0,0.0)
    V22 = [Point{3, Float64}(Q*v) for v ∈ V22] 
    V21 = V22 .+ Point{3, Float64}([0.0,-boxDim[2]/2,0.0])
    V22 = V22 .+ Point{3, Float64}([0.0, boxDim[2]/2,0.0])
    F21 = invert_faces(F22)

    Q = RotXYZ(0.0,0.5*pi,0.0)
    V32 = [Point{3, Float64}(Q*v) for v ∈ V32] 
    V31 = V32 .+ Point{3, Float64}([-boxDim[1]/2,0.0,0.0])
    V32 = V32 .+ Point{3, Float64}([ boxDim[1]/2,0.0,0.0])
    F31 = invert_faces(F32)

    F,V,C = joingeom(F11,V11,F12,V12,F21,V21,F22,V22,F31,V31,F32,V32)
    F,V,_,_ = mergevertices(F,V)
    return F,V,C
end

"""
    tetbox(boxDim,pointSpacing; stringOpt = "paAqYQ",region_vol=nothing)

Creates tetrahedral box mesh

# Description 
This function uses the dimensions defined in `boxDim`, and the desired point
spacing defined by `pointSpacing`, to create a tetrahedral mesh for a box. 
The output consists of the elements `E`, the vertices `V`, the boundary faces 
`Fb`, and the boundary face labelling `Cb` (for the 6 sides of the box). 
"""
function tetbox(boxDim,pointSpacing; stringOpt = "paAqYQ",region_vol=nothing)
    # Created triangulated box
    Fb,Vb,Cb = tribox(boxDim,pointSpacing)

    # Create tetrahedral mesh of the box
    E, V, _, Fb, Cb = tetgenmesh(Fb, Vb; facetmarkerlist=Cb, V_regions=nothing,region_vol=region_vol,V_holes=nothing,stringOpt=stringOpt)
    
    return E, V, Fb, Cb
end

"""
    pad3(A::Array{T,3}; padAmount = 1, padValue = T(0.0)) where T<:Real

Pads 3D array 

# Description 
This function pads the 3D input array `A` by the amount `padAmount` and with the
value `padValue`. The output is an array that is 2*padAmount larger in size 
direction. 
"""
function pad3(A::Array{T,3}; padAmount = 1, padValue = T(0.0)) where T<:Real
    siz = size(A) # Get size of A 

    # Create padded image 
    B = Array{T,3}(undef,siz .+ (2*padAmount)) 
    B[1:padAmount,:,:]   .= padValue
    B[:,1:padAmount,:]   .= padValue
    B[:,:,1:padAmount]   .= padValue
    B[end-padAmount+1:end,:,:] .= padValue
    B[:,end-padAmount+1:end,:] .= padValue
    B[:,:,end-padAmount+1:end] .= padValue

    # Set centre of B to A
    @inbounds for i in 1:1:siz[1]
        @inbounds for j in 1:1:siz[2]
            @inbounds for k in 1:1:siz[3]
                B[i+padAmount,j+padAmount,k+padAmount] = A[i,j,k]                
            end
        end
    end
    return B
end

"""
    getisosurface(A; level=0.0, cap=false, padValue=nothing, x::Union{AbstractVector{T},Nothing}=nothing, y::Union{AbstractVector{T},Nothing}, z::Union{AbstractVector{T},Nothing}) where T<:Real  

Constructs isosurface geometry 

# Description 
This function creates the triangular faces `F` and vertices `V` for an 
isosurface in the 3D image `A` of the level specified by `level`. 
"""
function getisosurface(A; level=0.0, cap=false, padValue=nothing, x::Union{AbstractVector{T},Nothing}=nothing, y::Union{AbstractVector{T},Nothing}=x, z::Union{AbstractVector{T},Nothing}=x) where T<:Real  
    if cap                
        # Get/determine padValue  
        if isnothing(padValue)
            if isapprox(level,0.0,atol=1e-8)
                padValue=1e10
            else
                padValue = level + 1*10^(round(Int,log10(abs(level)))+10)        
                if level<0                  
                    padValue *= 1.0
                end
            end
        end               
        A = pad3(A; padAmount = 1, padValue = padValue) # Pad input array
        if isnothing(x) || isnothing(y) || isnothing(z)            
            mc = MarchingCubes.MC(A)
        else
            # Extend coordinates to conform to padded array (assumes evenly spaced coordinates in each direction)
            s = x[2]-x[1]
            x = range(x[1]-s, x[end]+s, length(x)+2)            
            s = y[2]-y[1]
            y = range(y[1]-s, y[end]+s, length(y)+2)
            s = z[2]-z[1]
            z = range(z[1]-s, z[end]+s, length(z)+2)                        
            mc = MarchingCubes.MC(A; x=x, y=y, z=z)            
        end
    elseif !cap
        if isnothing(x) || isnothing(y) || isnothing(z)            
            mc = MarchingCubes.MC(A,Int)
        else
            mc = MarchingCubes.MC(A,Int; x=x, y=y, z=z)
        end
    end    
    MarchingCubes.march(mc,level)
    F = [TriangleFace{Int64}(f) for f in mc.triangles]
    V = [Point{3,Float64}(p) for p in mc.vertices]
    return F,V
end

"""
    randangle(siz::Union{Int,Tuple{Vararg{Int, N}}, Vector{Int}} = 1) where N

Returns random angles

# Description 
This function returns a random angle or array of random angles of the size 
`siz`. The angles are in radians and values lie between -pi and pi. 
"""
function randangle(siz::Union{Int,Tuple{Vararg{Int, N}}, Vector{Int}} = 1) where N
    if siz == 1
        return rand()*pi*rand((-1,1))
    else
        A = Array{Float64}(undef,siz)    
        for i in eachindex(A)
            A[i] = pi * rand() * rand((-1.0,1.0))
        end
        return A
    end    
end

"""
    stepfunc(type)

Returns a step function

# Description 
This function returns a step function (such as smoothstep functions [1]) to move
from a level `a` to `b` using the  function type specified by `type`. I.e. 
`f = stepfunc(type)` can be used as: `y = f(a,b,t)`
The functions are constrained such that the output is `a` if `t<=0.0` and the
output is `b` when `t>=1.0`. Each function uses the following definition: 
`(1.0 - f(t)) * a + f(t) * b`
Where `f(t)` depends on the function type requested. The following types are
supported: 
    :linear, this is a simple linear mapping (lerp) from `a` to `b`
    :Perlin, the Perlin smooth step function 6t⁵-15t⁴+10t³
    :smoothstep, 6t²-2t³
    :cosine, 2-cos(tπ)/2
The default is :linear

# References
1. https://en.wikipedia.org/wiki/Smoothstep
"""
function stepfunc(type::Symbol=:linear)      
    # Create inline function variations 
    f = t -> t    
    if type == :Perlin
        f = t -> t^3 * (t * (6.0 * t - 15.0) + 10.0) # 6t⁵-15t⁴+10t³ 
    elseif type == :smoothstep
        f = t -> 3.0 * t^2 - 2.0 * t^3 
    elseif type == :cosine
        f = t -> 0.5 - 0.5 * cos(t*pi) 
    elseif type == :linear
        # Nothing, already default
        # f = t -> t 
    else        
        error("InvalidParameter: Invalid type $type, valid options are :linear, :Perlin, :cosine, and :smoothstep")
    end            
    # Return function with constraints
    return function (a, b, t)
        if t <= 0.0
            return 0.0
        elseif t >= 1.0
            return 1.0
        else            
            return (1.0 - f(t)) * a + f(t) * b
        end
    end
end

"""
    perlin_noise(size_grid, sampleFactor, type=:Perlin)    

Returns Perlin noise array

# Description 
This function returns a 2D image containing Perlin noise [1]. The grid size is 
defined by `size_grid`. The `sampleFactor` defines the number of pixels to use 
for each grid cell. The output is a Matrix{Float64} with the size
`(size_grid .- 1) .* sampleFactor`. The `type` parameter dictates the type of 
"fade" function to use (see also: `stepfunc`), and the default is :Perlin. 

# References
1. https://en.wikipedia.org/wiki/Perlin_noise
"""
function perlin_noise(size_grid, sampleFactor; type=:Perlin)        
    pixelSize = 1/sampleFactor # Pixel size assuming grid has unit steps

    # Create grid vectors 
    A = randangle(size_grid) # Random angles    
    Ux = cos.(A) # Unit vector x-component
    Uy = sin.(A) # Unit vector y-component

    # Define "fade"/smoothstep function for interpolation 
    fade = stepfunc(type)

    # Initialise image
    size_image = (size_grid .- 1) .* sampleFactor # image size
    M = Matrix{Float64}(undef,size_image) # Start as undef Float64

    # Pre-compute grid cell quantities
    xy = range(0+pixelSize/2,1-pixelSize/2,sampleFactor) # x or y coordinates within a grid cell

    xc = [0,1,1,0]
    yc = [0,0,1,1]
    @inbounds for ip in 1:sampleFactor # For each cell row
        @inbounds for jp in 1:sampleFactor # For each cell column
            @inbounds for ig in 1:size_grid[1]-1 # For each grid row    
                @inbounds for jg in 1:size_grid[2]-1 # For each grid column
                    i = (ig-1)*sampleFactor + ip # Pixel row index
                    j = (jg-1)*sampleFactor + jp # Pixel column index
                    
                    # Current pixel cell coordinates
                    px = xy[jp]
                    py = xy[ip]
                            
                    # Offset vector components
                    xc1 = px # -xc[1] Offset vector 1 x
                    xc2 = px-xc[2] # Offset vector 2 x
                    xc3 = px-xc[3] # Offset vector 3 x
                    xc4 = px-xc[4] # Offset vector 4 x

                    yc1 = py # -yc[2] Offset vector 1 y                
                    yc2 = py-yc[2] # Offset vector 2 y               
                    yc3 = py-yc[3] # Offset vector 3 y                
                    yc4 = py-yc[4] # Offset vector 4 y

                    u1x = Ux[ig  , jg]
                    u2x = Ux[ig  , jg+1]
                    u3x = Ux[ig+1, jg+1]
                    u4x = Ux[ig+1, jg]

                    u1y = Uy[ig  ,jg]
                    u2y = Uy[ig  ,jg+1]
                    u3y = Uy[ig+1,jg+1]
                    u4y = Uy[ig+1,jg]

                    d1 = xc1.*u1x + yc1.*u1y
                    d2 = xc2.*u2x + yc2.*u2y
                    d3 = xc3.*u3x + yc3.*u3y
                    d4 = xc4.*u4x + yc4.*u4y
            
                    # Interpolation using fade function 
                    d12 = fade(d1,   d2, px)
                    d34 = fade(d4,   d3, px)
                    d   = fade(d12, d34, py)

                    M[i,j] = d                
                end
            end
        end
    end
    return M
end

"""
    removepoints(V,indRemove)

Removes listed points

# Description 
This function removes the points `indRemove` from the point vector `V`. The 
output contains a vector lacking these points, as well as `indFix`, a vector of
indices that can be used to update indices from before the point removal to 
equivalent indices after point removal. 
"""
function removepoints(V,indRemove)
    n = length(V)        
    indFix = zeros(Int,n)
    indKeep = setdiff(1:n,indRemove)
    indFix[indKeep] .= 1:length(indKeep)
    V = V[indKeep]    
    return V,indFix
end

"""
    inpolygon(p::Point{M,T}, V::Vector{Point{N,T}}; atol = sqrt(eps(T))) where N where M where T<:Real    
    inpolygon(P::Vector{Point{M,T}}, V::Vector{Point{N,T}}; atol = sqrt(eps(T))) where N where M where T<:Real

Finds points in/on/out of polygon

# Description 
This function is an improved version of the algorithm proposed by Hao et al. [1]
for determining if the input points `P` are inside the polygon defined by the 
input `V`. For a single point the output is a single integer which is: 
* -1 if the point is computed to be outside of the polygon
*  0 if the point is computed to be on the boundary of the polygon
*  1 if the point is computed to be inside the boundary
If the input is instead a vector of points then the output consists of a 
corresponding vector of such integers. 
This implementation differs from [1] in terms of the use of `atol` for 
approximate equivalance. This helps in more robust labelling of "on polygon" 
points.  

# References
1. https://doi.org/10.3390/sym10100477
"""
function inpolygon(p::Point{M,T}, V::Vector{Point{N,T}}; atol = sqrt(eps(T)), in_flag::TP=1, on_flag::TP=0, out_flag::TP=-1) where {TP} where N where M where T<:Real    
    k = 0
    xₚ = p[1]
    yₚ = p[2]    
    j = length(V)
    Ø = zero(T)
    for i in 1:1:length(V) # Looping over edges
        Δyᵢ = V[i][2] - yₚ # Point i y distance from p
        Δyⱼ = V[j][2] - yₚ # Point j y distance from p 
        isneg_Δyᵢ = Δyᵢ <= -atol # Point i significantly under p
        isneg_Δyⱼ = Δyⱼ <= -atol # Point j significantly under p       
        ispos_Δyᵢ = Δyᵢ >   atol # Point i significantly above p
        ispos_Δyⱼ = Δyⱼ >   atol # Point j significantly above p                       
        if ((isneg_Δyᵢ && isneg_Δyⱼ) || (ispos_Δyᵢ && ispos_Δyⱼ)) # Case 11 or 26 
            # edge is fully above or fully below
            j=i 
            continue
        end
        Δxᵢ = V[i][1] - xₚ # Point i x distance from p
        Δxⱼ = V[j][1] - xₚ # Point j x distance from p       
        if (isapprox(Δxᵢ,Ø,atol=atol) && isapprox(Δyᵢ,Ø,atol=atol)) || (isapprox(Δxⱼ,Ø,atol=atol) && isapprox(Δyⱼ,Ø,atol=atol))
            return on_flag
        end

        a = (Δxᵢ * Δyⱼ) - (Δxⱼ * Δyᵢ) # Edge cross product -> paralellogram area
        if ispos_Δyⱼ && isneg_Δyᵢ # Case 3, 9, 16, 21, 13, 24            
            if a > atol # Case 3 or 9
                k += 1
            elseif isapprox(a,Ø,atol=atol) # Case 16 or 21
                return on_flag 
            # else case 13 or 24
            end
        elseif ispos_Δyᵢ && isneg_Δyⱼ # Case 4, 10, 19, 20, 12, or 25            
            if a < -atol # Case 4 or 10 
                k += 1
            elseif isapprox(a,Ø,atol=atol) # Case 19 or 20
                return on_flag 
            # else case 12 or 25
            end
        elseif (isapprox(Δyⱼ,Ø,atol=atol) && Δyᵢ < -atol) # Case 7, 14, or 17            
            if a > atol #&& !((Δxⱼ <= -atol && Δxᵢ >= atol) || (Δxᵢ <= -atol && Δxⱼ >= atol))
                k += 1
            # elseif isapprox(a,Ø,atol=atol) # Case 17 -> only for point on vertex, handled already
            #     return on_flag             
            end   
        elseif (isapprox(Δyᵢ,Ø,atol=atol) && Δyⱼ < -atol) # Case 8, 15, or 18                           
            if a < -atol #&& !((Δxⱼ <= -atol && Δxᵢ >= atol) || (Δxᵢ <= -atol && Δxⱼ >= atol))
                k += 1
            # elseif isapprox(a,Ø,atol=atol) # Case 18 -> only for point on vertex, handled already
            #     return on_flag             
            end
        elseif isapprox(Δyᵢ,Ø,atol=atol) && isapprox(Δyⱼ,Ø,atol=atol) # Case 1, 2, 5, 6, 22, or 23                                       
            if (Δxⱼ <= -atol && Δxᵢ >= atol) # Case 1
                return on_flag
            elseif (Δxᵢ <= -atol && Δxⱼ >= atol) # Case 2
                return on_flag
                # else case 5, 6, 22, 23
            end
        end        
        j = i
    end
    iszero(k % 2) && return out_flag
    return in_flag
end

function inpolygon(P::Vector{Point{M,T}}, V::Vector{Point{N,T}}; atol = sqrt(eps(T)), in_flag::TP=1, on_flag::TP=0, out_flag::TP=-1) where {TP} where N where M where T<:Real
    return [inpolygon(p,V; atol = atol, in_flag=in_flag, on_flag=on_flag, out_flag=out_flag) for p in P]
end