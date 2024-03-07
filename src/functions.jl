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

This function simply returns the string for the Comodo path. This is helpful for instance to load items from the `assets`` folder. 
"""
function comododir()
    joinpath(@__DIR__, "..")
end


"""
    slidercontrol(hSlider,ax)    

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
    on(events(ax).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat # Pressed or held for instance
            if event.key == Keyboard.up                                  
                sliderValue = hSlider.value.val              
                if sliderValue ==hSlider.range.val[end]                
                    set_close_to!(hSlider, hSlider.range.val[1])
                else                
                    set_close_to!(hSlider, hSlider.value.val+1)
                end            
            elseif event.key == Keyboard.down
                sliderValue = hSlider.value.val            
                if sliderValue ==hSlider.range.val[1]
                    set_close_to!(hSlider, hSlider.range.val[end])
                else                
                    set_close_to!(hSlider, hSlider.value.val-1)
                end                        
            end
            if event.key == Keyboard.right                                              
                set_close_to!(hSlider, hSlider.value.val+1)            
            elseif event.key == Keyboard.left            
                set_close_to!(hSlider, hSlider.value.val-1)
            end
        end
    end    
end

"""
    elements2indices(F)

# Description
    
This function obtains the unique set of indices for the vertices (nodes) 
used by the the simplices defined by `F`. The vector `F` may contain any type of 
simplices. For instance the elements in `F` may be of the type 
`GeometryBasics.TriangleFace` or `GeometryBasics.QuadFace` (or any other) for 
surface mesh data. However volumetric elements of any type are permitted. In
essence this function simply returns `unique(reduce(vcat,F))`.
Hence any suitable vector containing vectors of numbers permitted by 
`reduce(vcat,F)` is supported. 
"""
function elements2indices(F)
    return unique(reduce(vcat,F))
end

"""
    gridpoints(x::Vector{T}, y=x, z=x) where T<:Real

# Description

The `gridpoints` function returns a vector of 3D points which span a grid in 3D 
space. Points are defined as per the input ranges or range vectors. The output 
point vector contains elements of the type `GeometryBasics.Point3`. 
"""
function gridpoints(x::Union{Vector{T}, AbstractRange{T}}, y=x, z=x) where T<:Real
    reshape([GeometryBasics.Point{3, T}(x, y, z) for z in z, y in y, x in x], 
                             length(x)*length(y)*length(z))
end  

"""
    interp_biharmonic_spline(x,y,xi; extrapolate_method=:linear,pad_data=:linear)

# Description

This function uses biharmonic spline interpolation. The input is assumed to 
represent ordered data representing a curve. 

# References
[David T. Sandwell, Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data, Geophysical Research Letters, 2, 139-142, 1987. doi: 10.1029/GL014i002p00139](https://doi.org/10.1029/GL014i002p00139)
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
        error("Invalid pad_data method provided, valued options are :linear, :constant, and :none")
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
        error("Invalid extrapolate_method method provided, valued options are :linear, :constant, and :biharmonic")
    end

    return yi
end


"""
    interp_biharmonic(x,y,xi)

# Description 

This function uses biharmonic interpolation. The input `x` should define a 
vector consisting of m points which are n-dimensional, and the input `y` should
be a vector consisting of m scalar data values. 

# References 

[David T. Sandwell, Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter data, Geophysical Research Letters, 2, 139-142, 1987. doi: 10.1029/GL014i002p00139](https://doi.org/10.1029/GL014i002p00139)
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

# Description 

This function returns `n` points for an m-th order Bézier spline, based on the 
m control points contained in the input vector `P`. This function supports point
vectors with elements of the type `AbstractPoint{3}` (e.g.
`GeometryBasics.Point{3, Float64}`) or `Vector{Float64}`.
"""
function nbezier(P::Vector{T},n::Integer) where T<:Union{AbstractPoint{3},Vector{Float64}}
    if n<2
        error("The vale of n is too low. Request at least two data points")
    end
    t = range(0,1,n) # t range
    N = length(P) # Number of control points 
    nn = 0:1:N-1
    nnr = N-1:-1:0
    f = factorial.(nn) 
    s = factorial(N-1)./( f.*reverse(f) ); # Sigma
    
    if T<:AbstractPoint{3}
        V =  [T(0.0,0.0,0.0) for _ ∈ 1:n]
    else
        V =  fill(zeros(Float64,3),n)
    end    
    
    for i ∈ 1:n
        b = s.* ((1.0-t[i]).^nnr) .* (t[i].^nn)
        for j = 1:1:N
            V[i] += P[j].*b[j]
        end
    end 
    return V
end

function lerp(x,y,xi)    
    if length(xi)>1
        yi=zeros(eltype(xi),length(xi))
        for q ∈ eachindex(xi)
            yi[q] = lerp_(x,y,xi[q])
        end
    else
        yi = lerp_(x,y,xi)
    end 
    return yi
end

function lerp_(x,y,xi)
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

function dist(V1,V2)
    D = Matrix{Float64}(undef,length(V1),length(V2))   
    for i ∈ eachindex(V1)
        for j ∈ eachindex(V2)          
            D[i,j] = euclidean(V1[i],V2[j]) # norm(V1[i]-V2[j])       
        end
    end
    return D
end

function dist(V1::Vector{T},v2::T) where T <: AbstractVector
    D = Matrix{Float64}(undef,length(V1),1)   
    for i ∈ eachindex(V1)        
        D[i,1] = euclidean(V1[i],v2)
    end
    return D
end

function dist(v1::T,V2::Vector{T}) where T <: AbstractVector
    D = Matrix{Float64}(undef,1,length(V2))   
    for j ∈ eachindex(V2)        
        D[1,j] = euclidean(v1,V2[j])
    end
    return D
end

function mindist(V1,V2; getIndex=false, skipSelf = false )
    D = Vector{Float64}(undef,length(V1))
    d = Vector{Float64}(undef,length(V2))
    if getIndex
        I = Vector{Int64}(undef,length(V1))
    end
    for i ∈ eachindex(V1)
        for j ∈ eachindex(V2)
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

function unique_dict_index(X::AbstractVector{T}; sort_entries=false) where T<:Real
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Nothing}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    for i ∈ eachindex(X)        
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

function unique_dict_index_inverse(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    indInverse = Vector{Int64}(undef,length(X)) 
    j=0
    for i ∈ eachindex(X)         
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

function unique_dict_index_count(X::AbstractVector{T}; sort_entries=false) where T<:Real
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    c =  Vector{Int64}() 

    j=0
    for i ∈ eachindex(X)      
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

function unique_dict_index_inverse_count(X::AbstractVector{T}; sort_entries=false) where T <: Real
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    indInverse = Vector{Int64}(undef,length(X)) 
    c =  Vector{Int64}() 

    j=0
    for i ∈ eachindex(X)      
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

function unique_dict_count(X::AbstractVector{T}; sort_entries=false) where T <: Real 
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    c = Vector{Int64}()
    j = 0
    for i ∈ eachindex(X)      
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

function unique_dict_inverse(X::AbstractVector{T}; sort_entries=false) where T <: Real 
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    # T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indInverse = Vector{Int64}(undef,length(X)) 

    j=0
    for i ∈ eachindex(X)      
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
    ind2sub(siz,ind)

# Description
Converts the linear indices in `ind`, for a matrix/array with size `siz`, to the 
equivalent subscript indices.  
"""
function ind2sub(siz,ind)
    
    numDim = length(siz)
    k = cumprod([siz[i] for i ∈ eachindex(siz)],dims=1)
    m = prod(siz)
    
    if any(ind.>m) || any(ind.<1)
        error("Encountered index value of out of valid range 1:$m")
    end

    if length(ind)>1
        A = [ind2sub_(ind_i,numDim,k) for ind_i ∈ ind]
    else
        A = ind2sub_(ind,numDim,k)      
    end

    return A
end

# ind2sub helper function to parse just a single linear index and produce a single subscript index set 
function ind2sub_(ind,numDim,k)

    a = Vector{Int64}(undef,numDim) # Initialise a

    for q ∈ numDim:-1:1   # For all dimensions     
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

# Description

Converts the subscript indices in `A`, for a matrix/array with size `siz`, to 
the equivalent linear indices.  
"""
function sub2ind(siz,A)
    
    numDim = length(siz)
    k = cumprod([siz[i] for i ∈ eachindex(siz)],dims=1)

    ind = Vector{Int64}(undef,length(A))
    for i ∈ eachindex(A)        
        ind[i] = round(Int64,A[i][1]);
        for q=2:1:numDim
            ind[i] += (A[i][q].-1).*k[q-1]
        end        
    end

    return ind
end

function sub2ind(siz,A::Vector{Int64})
    
    numDim = length(siz)
    k = cumprod([siz[i] for i ∈ eachindex(siz)],dims=1)
    
    ind = A[1];
    for q=2:1:numDim
        ind += (A[q].-1).*k[q-1]
    end        

    return ind
end

# Function to obtain the edges of a face based mesh
function meshedges(S; unique_only=false)
    # Get number of nodes per simplex, use first for now (better to get this from some property)
    s1 = S[1] # First simplex
    m = length(s1) #Number of nodes per simplex
    
    E = GeometryBasics.LineFace{Int}[]    
    for j1 ∈ 1:1:m # Loop over each node/point for the current simplex           
        if j1<m
            j2 = j1+1
        else
            j2 = 1
        end            
        for s ∈ S # Loop over each simplex        
            push!(E,(s[j1],s[j2]))            
        end 
    end

    # Remove doubles e.g. 1-2 seen as same as 2-1
    if unique_only
        E = gunique(E; sort_entries=true);
    end

    return E
end

function meshedges(M::GeometryBasics.Mesh; unique_only=false)   
    return meshedges(faces(M); unique_only=unique_only)
end

# Function to create the surface geometry data for an icosahedron
function icosahedron(r=1.0)

    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(ϕ + 2.0)
    t = ϕ*s

    # Create vertices
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,12)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s,  -t)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s,   t)
    V[3 ] = GeometryBasics.Point{3, Float64}( 0.0,    s,   t)
    V[4 ] = GeometryBasics.Point{3, Float64}( 0.0,    s,  -t)
    V[5 ] = GeometryBasics.Point{3, Float64}(  -s,   -t, 0.0)
    V[6 ] = GeometryBasics.Point{3, Float64}(  -s,    t, 0.0)
    V[7 ] = GeometryBasics.Point{3, Float64}(   s,    t, 0.0)
    V[8 ] = GeometryBasics.Point{3, Float64}(   s,   -t, 0.0)
    V[9 ] = GeometryBasics.Point{3, Float64}(  -t,  0.0,  -s)
    V[10] = GeometryBasics.Point{3, Float64}(   t,  0.0,  -s)
    V[11] = GeometryBasics.Point{3, Float64}(   t,  0.0,   s)
    V[12] = GeometryBasics.Point{3, Float64}(  -t,  0.0,   s)

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

# Function to create the surface geometry data for a octahedron
function octahedron(r=1.0)

    s = r/sqrt(2.0)

    # Create vertices    
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,6)
    V[1 ] = GeometryBasics.Point{3, Float64}(   -s,    -s,  0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}(    s,    -s,  0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(    s,     s,  0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}(   -s,     s,  0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}(  0.0,   0.0,   -r)
    V[6 ] = GeometryBasics.Point{3, Float64}(  0.0,   0.0,    r)
    
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

# Function to create the surface geometry data for a dodecahedron
function dodecahedron(r=1.0)

    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(3.0)
    t = ϕ*s    
    w = (ϕ-1.0)*s

    # Create vertices    
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,20)
    V[1 ] = GeometryBasics.Point{3, Float64}(   s,   s,   s)
    V[2 ] = GeometryBasics.Point{3, Float64}(   w, 0.0,   t)
    V[3 ] = GeometryBasics.Point{3, Float64}(  -t,  -w, 0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}(   t,   w, 0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}(  -s,   s,  -s)
    V[6 ] = GeometryBasics.Point{3, Float64}( 0.0,  -t,  -w)
    V[7 ] = GeometryBasics.Point{3, Float64}(  -t,   w, 0.0)
    V[8 ] = GeometryBasics.Point{3, Float64}(   s,  -s,   s)
    V[9 ] = GeometryBasics.Point{3, Float64}(  -s,   s,   s)
    V[10] = GeometryBasics.Point{3, Float64}(  -s,  -s,   s)
    V[11] = GeometryBasics.Point{3, Float64}(   s,  -s,  -s)
    V[12] = GeometryBasics.Point{3, Float64}(   w, 0.0,  -t)
    V[13] = GeometryBasics.Point{3, Float64}(  -s,  -s,  -s)
    V[14] = GeometryBasics.Point{3, Float64}( 0.0,  -t,   w)
    V[15] = GeometryBasics.Point{3, Float64}( 0.0,   t,  -w)
    V[16] = GeometryBasics.Point{3, Float64}(  -w, 0.0,   t)
    V[17] = GeometryBasics.Point{3, Float64}(   t,  -w, 0.0)
    V[18] = GeometryBasics.Point{3, Float64}(  -w, 0.0,  -t)
    V[19] = GeometryBasics.Point{3, Float64}(   s,   s,  -s)
    V[20] = GeometryBasics.Point{3, Float64}( 0.0,   t,   w)

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

# Function to create the surface geometry data for a cube
function cube(r=1.0)

    # Create vertices       
    s = r/sqrt(3.0)

    V=Vector{GeometryBasics.Point{3, Float64}}(undef,8)
    V[1 ] = GeometryBasics.Point{3, Float64}( -s, -s, -s)
    V[2 ] = GeometryBasics.Point{3, Float64}( -s,  s, -s)
    V[3 ] = GeometryBasics.Point{3, Float64}(  s,  s, -s)
    V[4 ] = GeometryBasics.Point{3, Float64}(  s, -s, -s)
    V[5 ] = GeometryBasics.Point{3, Float64}( -s, -s,  s)
    V[6 ] = GeometryBasics.Point{3, Float64}( -s,  s,  s)
    V[7 ] = GeometryBasics.Point{3, Float64}(  s,  s,  s)
    V[8 ] = GeometryBasics.Point{3, Float64}(  s, -s,  s)
        
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

# Function to create the surface geometry data for a tetrahedron
function tetrahedron(r=1.0)

    # Create vertices    
    a = r*sqrt(2.0)/sqrt(3.0)
    b = -r*sqrt(2.0)/3.0
    c = -r/3.0       

    V=Vector{GeometryBasics.Point{3, Float64}}(undef,4)
    V[1 ] = GeometryBasics.Point{3, Float64}(   -a,      b,   c)
    V[2 ] = GeometryBasics.Point{3, Float64}(    a,      b,   c)    
    V[3 ] = GeometryBasics.Point{3, Float64}(  0.0,    0.0,   r)
    V[4 ] = GeometryBasics.Point{3, Float64}(  0.0, -2.0*b,   c)  

    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,4)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    F[2 ] = TriangleFace{Int64}(4,2,1)
    F[3 ] = TriangleFace{Int64}(4,3,2)
    F[4 ] = TriangleFace{Int64}(4,1,3)
    
    return GeometryBasics.Mesh(V,F)
end

function togeometrybasics_faces(FM)
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)
    n = length(FM)
    m = length(FM[1])
    if m == 3 # Triangles
        F = [TriangleFace{Int64}(FM[q]) for q ∈ eachindex(FM)]        
    elseif m ==4 # Quads
        F = [QuadFace{Int64}(FM[q]) for q ∈ eachindex(FM)]        
    else # Other mesh type        
        F = [NgonFace{Int64}(FM[q]) for q ∈ eachindex(FM)]        
    end
    return F
end

function togeometrybasics_faces(FM::Matrix{Int64})
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)
    n, m = size(FM)
    if m == 3 # Triangles
        F = Vector{TriangleFace{Int64}}(undef, n)
        for q ∈ 1:1:n            
            F[q] = TriangleFace{Int64}(FM[q,:])
        end
    elseif m == 4 # Quads
        F = Vector{QuadFace{Int64}}(undef, n)
        for q ∈ 1:1:n            
            F[q] = QuadFace{Int64}(FM[q,:])
        end
    else # Other mesh type
        F = Vector{m,NgonFace{Int64}}(undef, n)        
        for q ∈ 1:1:n            
            F[q] = NgonFace{m,Int64}(FM[q,:])
        end
    end
    return F
end

function togeometrybasics_points(VM)
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    n = length(VM)
    V=Vector{GeometryBasics.Point{3, Float64}}(undef, n)
    for q ∈ 1:1:n
        V[q] = GeometryBasics.Point{3, Float64}(VM[q])
    end
    return V
end

function togeometrybasics_points(VM::Matrix{Float64})
    n = length(VM)
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef, n)
    for q ∈ 1:1:n
        V[q] = GeometryBasics.Point{3, Float64}(VM[q,:])
    end
    return V
end

function togeometrybasics_points(VM::Vector{Vector{Int64}})
    n = length(VM)
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef, n)
    for q ∈ 1:1:n
        V[q] = GeometryBasics.Point{3, Float64}(VM[q])
    end
    return V
end

function togeometrybasics_mesh(VM,FM)
    V=togeometrybasics_points(VM)
    F=togeometrybasics_faces(FM)
    return GeometryBasics.Mesh(V,F)
end


function edgecrossproduct(F,V)    
    C = Vector{GeometryBasics.Vec{3, Float64}}(undef,length(F)) # Allocate array cross-product vectors
    n =  length(F[1]) # Number of nodes per face    
    for q ∈ eachindex(F) # Loop over all faces
        c  = cross(V[F[q][n]],V[F[q][1]]) # Initialise as cross product of last and first vertex position vector
        for qe ∈ 1:1:n-1 # Loop from first to end-1            
            c  += cross(V[F[q][qe]],V[F[q][qe+1]]) # Add next edge contribution          
        end
        C[q] = c./2 # Length = face area, direction is along normal vector
    end
    return C
end

function edgecrossproduct(M::GeometryBasics.Mesh) 
    return edgecrossproduct(faces(M),coordinates(M))
end

function facenormal(M::GeometryBasics.Mesh) 
    return facenormal(faces(M),coordinates(M))
end

# Computes mesh face normals
function facenormal(F,V)
    C = edgecrossproduct(F,V)
    return C./norm.(C)
end

# Computes mesh face areas
function facearea(F,V)        
    return norm.(edgecrossproduct(F,V))
end

function facearea(M::GeometryBasics.Mesh)         
    return facearea(faces(M),coordinates(M))
end

function vertexnormal(M::GeometryBasics.Mesh; weighting=:area)     
    return normalizevector(vertexnormal(faces(M),coordinates(M); weighting=weighting))
end

function vertexnormal(F,V; weighting=:area)     
    return normalizevector(simplex2vertexdata(F,facenormal(F,V),V; weighting=weighting))
end

function edgelengths(F,V)    
    return [norm(V[e[1]]-V[e[2]]) for e ∈ meshedges(F; unique_only=true)]
end


"""
    platonicsolid(n,r=1.0)

# Description

Creates a GeometryBasics mesh description for a platonic solid of choice. The 
input `n` defines the choice.
1. tetrahedron
2. cube
3. octahedron
4. icosahedron
5. dodecahedron

The final input parameter `r` defines the radius of the platonic solid (the 
radius of the circumsphere to the vertices).

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

function unique_dict(X::AbstractVector{T}) where T <: Real    
    d = OrderedDict{T ,Int64}() # Use dict to keep track of used values    
    indUnique = Vector{Int64}()
    indReverse = Vector{Int64}(undef,length(X))
    # c = Vector{Int64}()
    j=0
    for i ∈ eachindex(X)        
        if !haskey(d, X[i])                                          
            j+=1
            d[X[i]] = j # reverse index in dict            
            push!(indUnique, i)                     
            indReverse[i] = j 
        else
            indReverse[i] = d[X[i]]            
        end        
    end
    return collect(keys(d)), indUnique, indReverse #, c
end

function unique_simplices(F,V=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,F)) 
    else
        n = length(V)
    end

    virtualFaceIndices = sub2ind(n.*ones(Int64,length(F[1])),sort.(F))    
    ~, ind1, ind2 = unique_dict(virtualFaceIndices) 

    return F[ind1], ind1, ind2
end

function subtri(F,V,n; method = :linear)
    
    if iszero(n)
        return F,V
    elseif isone(n)
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        Fm1 = [TriangleFace{Int64}(a.+length(V)) for a ∈ eachrow(reshape(indReverse,length(F),length(F[1])))] 
        Fm2 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm3 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm4 = Vector{TriangleFace{Int64}}(undef,length(Fm1))        
        for i ∈ eachindex(F)                        
            Fm2[i] = TriangleFace{Int64}([Fm1[i][1], Fm1[i][3], F[i][1]])
            Fm3[i] = TriangleFace{Int64}([Fm1[i][2], Fm1[i][1], F[i][2]])
            Fm4[i] = TriangleFace{Int64}([Fm1[i][3], Fm1[i][2], F[i][3]])
        end

        # Create combined face set
        Fn = [Fm1; Fm2; Fm3; Fm4]        

        con_E2F = con_edge_face(F,Eu,indReverse)

        # Create new vertices depending on method
        if method == :linear # Simple linear splitting
            # Create complete point set
            Vn = [V; simplexcenter(Eu,V)]  # Old and new mid-edge points          
        elseif method == :loop #Loop subdivision 
    
            # New mid-edge like vertices
            Vm = Vector{GeometryBasics.Point{3, Float64}}(undef,length(Eu)) 
            for q ∈ eachindex(Eu) # For each edge index                        
                F_touch = F[con_E2F[q]] # Faces sharing current edge, mostly 2 but 1 for a boundary edge
                indVerticesTouch = Vector{Int64}() 
                for f ∈ F_touch        
                    b = f.!=Eu[q][1] .&& f.!=Eu[q][2]      
                    if any(b)  
                        append!(indVerticesTouch,f[b])           
                    end
                end        
                Vm[q]=3/8 .*(V[Eu[q][1]] .+ V[Eu[q][2]])  .+ 1/8 .* (V[indVerticesTouch[1]] .+ V[indVerticesTouch[2]])
            end
    
            # Modified vertices for original vertices
            Vv = Vector{GeometryBasics.Point{3, Float64}}(undef,length(V))
            for q ∈ eachindex(V)            
                B_vert_face = [any(f.==q) for f in F]
                F_touch = F[B_vert_face] # Faces mostly 2 but 1 for a boundary edge
                indVerticesTouch = Vector{Int64}()
                for f ∈ F_touch                
                    indTouch = f[f.!=q]        
                    for i ∈ indTouch 
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
            error("Incorrect metod. Use: :linear or :loop")
        end

        return Fn,Vn    

    elseif n>1
        for _ =1:1:n
            F,V = subtri(F,V,1; method=method)
        end
        return F,V
    end
end

function subquad(F,V,n; method=:linear)
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
            Ve = Vector{GeometryBasics.Point{3, Float64}}(undef,length(Eu)) # Initialize edge points
            for q ∈ eachindex(Eu)                         
                Ve[q] = (mean(Vf[con_E2F[q]],dims=1)[1] .+ Ve_mid[q])./2.0
            end

            # Vertex points 
            Vv = Vector{GeometryBasics.Point{3, Float64}}(undef,length(V)) # Initialize vertex points
            for q ∈ eachindex(V) # Loop over all vertices
                indF = con_V2F[q]
                indE = con_V2E[q]
                N = length(indF) # Number of faces (or edges) touching this vertex                    
                Vv[q] = (mean(Vf[indF],dims=1)[1] .+ 2.0.*mean(Ve_mid[indE],dims=1)[1] .+ (N-3.0).*V[q])./N
            end

            Vn = [Vv;Ve;Vf] # Joined point set
        end

        # Define faces
        Fn = Vector{QuadFace{Int64}}(undef,length(F)*4)        
        nv = length(V)
        ne = length(Eu)
        for q ∈ eachindex(F)
            i = 1 + (q-1)*4
            for ii = 0:1:3
                Fn[i+ii] = QuadFace{Int64}([F[q][ii+1], con_F2E[q][ii+1]+nv, q+nv+ne, con_F2E[q][1+mod(3+ii,4)]+nv])                
            end            
        end
        return Fn,Vn
    elseif n>1
        for _ =1:1:n
            F,V = subquad(F,V,1;method=method)
        end

        return F,V
    end
end

# Create geodesic dome
function geosphere(n,r)
    M = platonicsolid(4,r)
    V = coordinates(M)
    F = faces(M)
    for _ = 1:n
        F,V = subtri(F,V,1)
        for q ∈ eachindex(V)
            v = V[q]
            rn = sqrt(sum(v.^2))
            V[q] = v .* (r/rn)
        end
    end
    return F,V
end

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

    for q ∈ 1:numElements
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

    V = convert(Vector{GeometryBasics.Point{3, Float64}},IJK_nodes)
    for q ∈ eachindex(V)
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

function con_face_edge(F,E_uni=missing,indReverse=missing)
    if ismissing(E_uni) | ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    return [Vector{Int64}(a) for a ∈ eachrow(reshape(indReverse,length(F),length(F[1])))] # [indReverse[[1,2,3].+ (i-1)*3] for i ∈ eachindex(F)]
end

function con_edge_face(F,E_uni=missing,indReverse=missing)
    if ismissing(E_uni) || ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    con_F2E = con_face_edge(F,E_uni,indReverse)
    
    con_E2F = [Vector{Int64}() for _ ∈ 1:1:length(E_uni)]
    for i_f ∈ eachindex(F)
        for i ∈ con_F2E[i_f]
            push!(con_E2F[i],i_f)
        end
    end 
    return con_E2F
end

function con_face_face(F,E_uni=missing,indReverse=missing,con_E2F=missing,con_F2E=missing)
    if ismissing(E_uni)| ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end    
    if ismissing(con_E2F) 
        con_E2F = con_edge_face(F,E_uni)
    end
    if ismissing(con_F2E)
        con_F2E = con_face_edge(F,E_uni,indReverse)                 
    end
    con_F2F = [Vector{Int64}() for _ ∈ 1:1:length(F)]
    for i_f ∈ eachindex(F)
        for i ∈ reduce(vcat,con_E2F[con_F2E[i_f]])    
            if i!=i_f     
                push!(con_F2F[i_f],i)
            end 
        end
    end
    return con_F2F
end

function con_face_face_v(F,V=missing,con_V2F=missing)
    if ismissing(con_V2F) 
        con_V2F = con_vertex_face(F,V)  # VERTEX-FACE connectivity
    end
    con_F2F = [Vector{Int64}() for _ ∈ 1:1:length(F)]
    for i_f ∈ eachindex(F)
        for i ∈ reduce(vcat,con_V2F[F[i_f]])    
            if i!=i_f     
                push!(con_F2F[i_f],i)
            end 
        end
    end
    return con_F2F
end

function con_vertex_simplex(F,V=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,F))
    else
        n = length(V)
    end
    con_V2F = [Vector{Int64}() for _ ∈ 1:1:n]
    for i_f ∈ eachindex(F)
        for i ∈ F[i_f]
            push!(con_V2F[i],i_f)
        end
    end
    return con_V2F
end

function con_vertex_face(F,V=missing)
    return con_vertex_simplex(F,V)
end

function con_vertex_edge(E,V=missing)
    return con_vertex_simplex(E,V)
end

function con_edge_edge(E_uni,con_V2E=missing)
    if ismissing(con_V2E)
        con_V2E = con_vertex_edge(E_uni) 
    end    
    con_E2E = [Vector{Int64}() for _ ∈ 1:1:length(E_uni)]
    for i_e ∈ eachindex(E_uni)
        for i ∈ reduce(vcat,con_V2E[E_uni[i_e]])    
            if i!=i_e     
                push!(con_E2E[i_e],i)
            end 
        end
    end
    return con_E2E
end

function con_vertex_vertex_f(F,V=missing,con_V2F=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,E))
    else 
        n = length(V)
    end
    
    if ismissing(con_V2F)
        con_V2F = con_vertex_face(F,V)
    end

    con_V2V = [Vector{Int64}() for _ ∈ 1:1:n]
    for i_v ∈ 1:1:n
        for i ∈ unique(reduce(vcat,F[con_V2F[i_v]]))
            if i_v!=i
                push!(con_V2V[i_v],i)
            end
        end
    end

    return con_V2V
end

function con_vertex_vertex(E,V=missing,con_V2E=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,E))
    else 
        n = length(V)
    end
    if ismissing(con_V2E)
        con_V2E = con_vertex_edge(E,V)
    end

    con_V2V = [Vector{Int64}() for _ ∈ 1:1:n]
    for i_v ∈ 1:1:n
        for i ∈ reduce(vcat,E[con_V2E[i_v]])
            if i_v!=i
                push!(con_V2V[i_v],i)
            end
        end
    end

    return con_V2V
end

function meshconnectivity(F,V) #conType = ["ev","ef","ef","fv","fe","ff","ve","vf","vv"]

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

function mergevertices(F,V; roundVertices = true, numDigitsMerge=missing)

    m = length(V)
    if roundVertices
        if ismissing(numDigitsMerge)
            E = meshedges(F)
            d = [sqrt( sum((V[e[1]] .- V[e[2]]).^2) ) for e ∈ E]
            pointSpacing = mean(d)
            m = round(Int64,log10(pointSpacing))
            numDigitsMerge = 6-m
        end

        # Create rounded coordinates to help obtain unique set
        VR = [round.(v,digits = numDigitsMerge) for v ∈ V]

        # Get unique indices and reverse for rounded vertices
        _,ind1,ind2 = gunique(VR; return_index=true, return_inverse=true,sort_entries=false)
        V = V[ind1] # The unique node set
    else
        V,ind1,ind2 = gunique(V; return_index=true, return_inverse=true,sort_entries=false)
    end

    if length(V) != m # If the length has changed
        # Correct indices for faces
        for q ∈ eachindex(F)
            F[q] = ind2[F[q]]
        end
    end

    return F,V,ind1,ind2
end

"""

    smoothmesh_laplacian(F,V,con_V2V=missing; n=1, λ=0.5)

# Description

This function implements Weighted Laplacian mesh smoothing. At each 
iteration, this method replaces each point by the weighted sum of the 
Laplacian mean for the point and the point itself. The weighting is 
controlled by the parameter λ which is in the range (0,1). If λ=0 no 
smoothing occurs. If λ=0 pure Laplacian mean based smoothing occurs. For 
intermediate values a linear blending between the two occurs.  
"""
function smoothmesh_laplacian(F,V,con_V2V=missing; n=1, λ=0.5)

    if λ!=1
        # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if missing
        if ismissing(con_V2V)
            E = meshedges(F)
            E_uni,_,_ = unique_simplices(E)
            con_V2V = con_vertex_vertex(E_uni)
        end
        for _ = 1:1:n
            Vs = deepcopy(V)
            for q ∈ eachindex(V)
                # Linear blend between original and pure Laplacian
                Vs[q] = (1.0-λ).*Vs[q] .+ λ*mean(V[con_V2V[q]])
            end
            V = Vs
        end
    end
    return V
end


"""
    smoothmesh_hc(F,V, con_V2V=missing; n=1, α=0.1, β=0.5, tolDist=missing)

# Description 

This function implements HC (Humphrey's Classes) smoothing. This method uses
Laplacian like smoothing but aims to compensate for the shrinkage/swelling 
seen with pure Laplacian smoothing. 

# Reference 

[Vollmer et al. Improved Laplacian Smoothing of Noisy Surface Meshes, 1999. doi: 10.1111/1467-8659.00334](https://doi.org/10.1111/1467-8659.00334)


"""
function smoothmesh_hc(F,V, con_V2V=missing; n=1, α=0.1, β=0.5, tolDist=missing)

    # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if missing
    if ismissing(con_V2V)
        E = meshedges(F)
        E_uni,_,_ = unique_simplices(E)
        con_V2V = con_vertex_vertex(E_uni)
    end
    P = deepcopy(V) # Copy original input points
    B = deepcopy(V) # Initialise B
    c = 0
    while c<n       
        Q = deepcopy(P) # Reset Q as P for this iteration
        for i ∈ eachindex(V)
            P[i] = mean(Q[con_V2V[i]]) # Laplacian 
            # Compute different vector between P and a point between original 
            # point and Q (which is P before laplacian)
            B[i] = P[i] .- (α.*V[i] .+ (1-α).*Q[i])
        end
        d = 0.0        
        for i ∈ eachindex(V)      
            # Push points back based on blending between pure difference vector
            # B and the Laplacian mean of these      
            P[i] = P[i] .- (β.*B[i] .+ (1-β).* mean(B[con_V2V[i]]))
        end   
        c+=1 
        if !ismissing(tolDist) # Include tolerance based termination
            d = 0.0
            for i ∈ eachindex(V)
                d+=sqrt(sum((P[i].-Q[i]).^2)) # Sum of distances
            end
            if d<tolDist # Sum of distance smaller than tolerance?
                break
            end            
        end
    end
    return P
end

function quadplate(plateDim,plateElem)
    num_x = plateElem[1]+1
    num_y = plateElem[2]+1
    V = Vector{GeometryBasics.Point{3, Float64}}()
    for y = range(-plateDim[2]/2,plateDim[2]/2,num_y)
        for x = range(-plateDim[1]/2,plateDim[1]/2,num_x)
            push!(V,GeometryBasics.Point{3, Float64}(x,y,0.0))
        end
    end

    F = Vector{QuadFace{Int64}}()
    ij2ind(i,j) = i + ((j-1)*num_x) # function to convert subscript to linear indices

    for i = 1:1:plateElem[1]
        for j = 1:1:plateElem[2]        
            push!(F,QuadFace{Int64}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)]))
        end
    end
    return F, V
end

function quadsphere(n,r)
    M = platonicsolid(2,r)
    F = faces(M)
    V = coordinates(M)
    if n > 0
        for q ∈ 1:1:n
            F,V = subquad(F,V,1)
            V = r .* (V ./ norm.(V))
        end
    end
    return F,V
end

function simplex2vertexdata(F,DF,V=missing; con_V2F=missing, weighting=:none)
    
    if ismissing(con_V2F)
        con_V2F = con_vertex_face(F,V)
    end
    if weighting==:area
        if ismissing(V)
           error("Vertices need to be provided for area based weighting.") 
        else
            A = facearea(F,V)
        end
    end    
    DV = (typeof(DF))(undef,length(con_V2F))
    T = eltype(DV)
    for q ∈ eachindex(DV)
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
    for q ∈ eachindex(F)
        DF[q] = mean(T,DV[F[q]]) # The mean of the vertex data for each entry in F
    end
    return DF
end

function simplexcenter(F,V)
    return vertex2simplexdata(F,V)
end

function normalizevector(A)
    if eltype(A) <: Number
        return A./norm(A)
    else
        return A./norm.(A)
    end
end

function circlepoints(r,n; dir=:acw)
    if dir==:acw
        return [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),0) for t ∈ range(0,2*π-(2*π)/n,n)]
    elseif dir==:cw
        return [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),0) for t ∈ range(0,(2*π)/n-2*π,n)]
    else
        error("Invalid dir specified, use :acw or :cw")
    end
end

function circlepoints(f::FunctionType,n; dir=:acw) where {FunctionType <: Function}
    if dir==:acw
        return [GeometryBasics.Point{3, Float64}(f(t)*cos(t),f(t)*sin(t),0) for t ∈ range(0,2*π-(2*π)/n,n)]
    elseif dir==:cw
        return [GeometryBasics.Point{3, Float64}(f(t)*cos(t),f(t)*sin(t),0) for t ∈ range(0,(2*π)/n-2*π,n)]
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
function loftlinear(V1,V2;num_steps=2,close_loop=true,face_type=:tri)

    num_loop = length(V1)
    T = eltype(V1)
    # Linearly blending points from first to last
    V = Vector{T}()
    for q ∈ range(0,num_steps,num_steps)
        λ = q/num_steps
        Vn = (1.0-λ).*V1 .+ λ.* V2 
        append!(V,Vn)
    end

    if face_type == :tri
        V0 = deepcopy(V)
        for qq ∈ 2:2:num_steps-1
            i = (1:1:num_loop) .+ (qq-1) *num_loop
            for q ∈ 1:1:num_loop    
                if q == num_loop       
                    V[i[q]] = 0.5 .* (V0[i[q]]+V0[i[1]]) 
                else
                    V[i[q]] = 0.5 .* (V0[i[q]]+V0[i[q]+1])
                end
            end
        end 
    end

    ij2ind(i,j) = i + ((j-1)*num_loop) # function to convert subscript to linear indices

    # Build faces
    if face_type == :quad    
        F = Vector{QuadFace{Int64}}()
        for i = 1:1:num_loop-1
            for j = 1:1:num_steps-1    
                push!(F,QuadFace{Int64}([ij2ind(i,j+1),ij2ind(i+1,j+1),ij2ind(i+1,j),ij2ind(i,j)  ]))
            end
        end

        # Add faces to close over shape if requested
        if close_loop
            for q ∈ 1:1:num_steps-1                
                push!(F,QuadFace{Int64}([ ij2ind(num_loop,q+1), ij2ind(1,q+1), ij2ind(1,q), ij2ind(num_loop,q) ])) 
            end
        end
    elseif face_type == :tri_slash 
        F = Vector{TriangleFace{Int64}}()
        for i = 1:1:num_loop-1
            for j = 1:1:num_steps-1    
                push!(F,TriangleFace{Int64}([ ij2ind(i+1,j+1), ij2ind(i+1,j), ij2ind(i,j)     ])) # 1 2 3
                push!(F,TriangleFace{Int64}([ ij2ind(i,j),     ij2ind(i,j+1), ij2ind(i+1,j+1) ])) # 3 4 1
            end
        end

        # Add faces to close over shape if requested
        if close_loop
            for q ∈ 1:1:num_steps-1
                push!(F,TriangleFace{Int64}([ ij2ind(1,q+1),      ij2ind(1,q),          ij2ind(num_loop,q)  ])) # 1 2 3
                push!(F,TriangleFace{Int64}([ ij2ind(num_loop,q), ij2ind(num_loop,q+1), ij2ind(1,q+1)       ])) # 3 4 1
            end
        end
    elseif face_type == :tri 
        F = Vector{TriangleFace{Int64}}()
        for i = 1:1:num_loop-1
            for j = 1:1:num_steps-1    
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
            for q ∈ 1:1:num_steps-1
                if iseven(q) 
                    push!(F,TriangleFace{Int64}([ ij2ind(num_loop,q), ij2ind(num_loop,q+1), ij2ind(1,q+1)      ])) 
                    push!(F,TriangleFace{Int64}([ ij2ind(1,q+1),      ij2ind(1,q),          ij2ind(num_loop,q) ])) 
                else
                    push!(F,TriangleFace{Int64}([ ij2ind(num_loop,q+1), ij2ind(1,q+1),      ij2ind(1,q)          ]))
                    push!(F,TriangleFace{Int64}([ ij2ind(1,q),          ij2ind(num_loop,q), ij2ind(num_loop,q+1) ]))   
                end
            end
        end
    else
        error("Invalid face_type specified, use :tri, :tri_slash, or :quad")
    end
    return F, V
end 

function dirplot(ax,V,U; color=:black,linewidth=3,scaleval=1.0,style=:from)
    E = [GeometryBasics.LineFace{Int}(i,i+length(V)) for i ∈ 1:1:length(V)]
    if style==:from
        P = vcat(V,V.+(scaleval.*U))
    elseif style==:to
        P = vcat(V.-(scaleval.*U),V)
    elseif style==:through
        UU = (scaleval.*U)/2
        P = vcat(V.-UU,V.+UU)
    else
        error("Invalid style specified, use :from, :to, or :through")
    end
    hp = wireframe!(ax,GeometryBasics.Mesh(P,E),linewidth=linewidth, transparency=false, color=color)
    return hp
end

function normalplot(ax,M; type_flag=:face, color=:black,linewidth=3,scaleval=missing)
    F = faces(M)
    V = coordinates(M)
    E = meshedges(F)
    if ismissing(scaleval)
        scaleval = mean([norm(V[e[1]]-V[e[2]])/2.0 for e ∈ E])
    end
    if type_flag == :face        
        N = facenormal(F,V)
        V = simplexcenter(F,V)        
    elseif type_flag == :vertex
        N = vertexnormal(F,V)          
    else
        error("Incorrect type_flag, use :face or :vertex")
    end 
    
    hp = dirplot(ax,V,N; color=color,linewidth=linewidth,scaleval=scaleval,style=:from)

    # hp = arrows!(ax,V,N.*d,color=color,quality=6)    
    return hp 
end

function wrapindex(i::UnitRange{Int64},n)
    return 1 .+ mod.(i .+ (n-1),n)
end

function wrapindex(i::Vector{Int64},n)
    return 1 .+ mod.(i .+ (n-1),n)
end

function wrapindex(i::Int64,n)
    return 1+mod(i+(n-1),n)
end

function edgeangles(F,V)
    # TO DO: Fix vector type for variable `a` below
    m = length(F[1])
    A = Vector{GeometryBasics.Vec{m, Float64}}()
    for f ∈ F
        a = Vector{Float64}(undef,m)
        for i ∈ 1:1:m                        
            ip1 = wrapindex(i+1,m)            
            ip2 = wrapindex(i+2,m)
            n1 = normalizevector(V[f[ip1]]-V[f[i]])
            n2 = normalizevector(V[f[ip2]]-V[f[ip1]])
            a[i] = acos( dot(n1,n2) )
        end
        push!(A,a)
    end
    return A
end

function quad2tri(F,V; convert_method = :angle)
    # Local functions for slash based conversion   
    forw_slash(f) = [[f[1],f[2],f[3]],[f[3],f[4],f[1]]] # Forward slash 
    back_slash(f) = [[f[1],f[2],f[4]],[f[2],f[3],f[4]]] # Back slash

    Ft = Vector{TriangleFace{Int64}}()#(undef,length(Fn1)*2)
    for f ∈ F        
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
            δaf = sum(abs.(reduce(vcat,af) .- (π/3)))
            δab = sum(abs.(reduce(vcat,ab) .- (π/3)))
            if δaf<δab
                ft = ff
            else
                ft = fb
            end    
        else
            error("Incorrect conver_method set, use :forward, :backward, or :angle")
        end
        push!(Ft,ft[1])
        push!(Ft,ft[2])
    end
    return Ft
end

function remove_unused_vertices(F,V)
    T = eltype(F)
    indUsed = elements2indices(F)
    Vc = V[indUsed] # Remove unused points    
    indFix = zeros(length(V))
    indFix[indUsed].=1:1:length(indUsed)
    Fc = [(T)(indFix[f]) for f ∈ F] # Fix indices in F 
    return Fc, Vc, indFix
end

function trisurfslice(F,V,n = (0.0,0.0,1.0), p = mean(V,dims=1); snapTolerance = 0, output_type=:full)

    intersectFunc(v1,v2,d,n) = v1 .- d/dot(n,v2.-v1) .* (v2.-v1)
    
    # Compute dot product with normal of slicing plane
    d = map(v-> dot(n,v.-p),V)
    if snapTolerance != 0.0
        d[abs.(d).<snapTolerance] .= 0
    end
    LV = d.<0
    
    Fn =  Vector{TriangleFace{Int64}}()
    Vn = deepcopy(V)
    D = Dict{Vector{Int64},Int64}() # For pointing from edge to intersection point index
    for f ∈ F
        lf = LV[f]
        
        if any(lf) # Some or all below
            if all(lf) # All below
                if output_type == :full || output_type == :below            
                    push!(Fn,f)
                end                 
            else # Some below -> cut
                nBelow = sum(lf) # Number of nodes below
                if isone(nBelow) # 2-above, 1 below 
                    indP = f[wrapindex(findfirst(lf) .+ (0:2),3)]
                    
                    e1 = sort(indP[[1,2]])
                    if ~haskey(D,e1)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[2]],d[indP[1]],n))
                        D[e1] = length(Vn)
                    end
                    
                    e2 = sort(indP[[1,3]])
                    if ~haskey(D,e2)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[3]],d[indP[1]],n))
                        D[e2] = length(Vn)
                    end

                    if output_type == :above || output_type == :full                        
                        push!(Fn,TriangleFace{Int64}(D[e1],indP[2],indP[3]))
                        push!(Fn,TriangleFace{Int64}(D[e1],indP[3],D[e2]))
                    elseif output_type == :below || output_type == :full
                        push!(Fn,TriangleFace{Int64}(indP[1],D[e1],D[e2]))
                    end

                else # 1-above, 2 below
                    indP = f[wrapindex(findfirst(.~lf) .+ (0:2),3)]

                    e1 = sort(indP[[1,2]])
                    if ~haskey(D,e1)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[2]],d[indP[1]],n))
                        D[e1] = length(Vn)
                    end                    

                    e2 = sort(indP[[1,3]])
                    if ~haskey(D,e2)
                        push!(Vn,intersectFunc(Vn[indP[1]],Vn[indP[3]],d[indP[1]],n))
                        D[e2] = length(Vn)
                    end
                    
                    if output_type == :below || output_type == :full                        
                        push!(Fn,TriangleFace{Int64}(D[e1],indP[2],indP[3]))
                        push!(Fn,TriangleFace{Int64}(D[e1],indP[3],D[e2]))
                    elseif output_type == :above || output_type == :full
                        push!(Fn,TriangleFace{Int64}(indP[1],D[e1],D[e2]))
                    end
                end
            end
        else # Not any below -> all above
            if output_type == :full || output_type == :above            
                push!(Fn,f)
            end    
        end
    end    
    return remove_unused_vertices(Fn,Vn)
end

function count_edge_face(F,E_uni=missing,indReverse=missing)
    if ismissing(E_uni) || ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    con_F2E = con_face_edge(F,E_uni,indReverse)
    
    C = zeros(Int64,length(E_uni))
    for i_f ∈ eachindex(F)
        for i ∈ con_F2E[i_f]
            C[i]+=1
        end
    end 
    return C
end

function boundaryedges(F)
    E = meshedges(F)
    Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
    count_E2F = count_edge_face(F,Eu,indReverse)
    return Eu[isone.(count_E2F)]
end

function edges2curve(Eb)
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
        elseif Eb[e_ind[2]][1]==ind[end] #Check if 1st point of 2nd edge equals end
            i = e_ind[2]
        end
    end
    return ind
end

function pointspacingmean(V)
    # Equivalent to:  mean(norm.(diff(V,dims=1)))
    p = 0.0
    n = length(V)
    for i ∈ 1:1:n-1
        p += norm(V[i].-V[i+1])/(n-1)
    end
    return p
end

function extrudecurve(V1,d; s=1, n=Point{3, Float64}(0.0,0.0,1.0),num_steps=missing,close_loop=false,face_type=:quad)
    if ismissing(num_steps)
        num_steps = ceil(Int64,d/pointspacingmean(V1))
        if face_type==:tri
            num_steps = num_steps + Int64(iseven(num_steps)) # Force uneven
        end
    end
    if isone(s) # Allong n from V1
        p = d.*n
    elseif isone(-s) # Against n from V1
        p = -d.*n
    elseif iszero(s) # Extrude both ways from V1
        p = d.*n
        V1 = [(eltype(V1))(v.-p./2) for v ∈ V1] #Shift V1 in negative direction
    end
    V2 = [(eltype(V1))(v.+p) for v ∈ V1]  
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
        error("Wrong `con_type`, use :v or :e")
    end

    if all(isempty.(con_F2F)) # Completely disconnected face set (e.g. raw STL import)
        C = collect(1:1:length(F))
    else
        C = fill(0,length(F))
        i = 1
        c = 1
        C[i] = c 
        seen = Set{Int64}(1)
        while length(seen)<length(F)
            np = length(seen)
            con_f2f = con_F2F[i]
            if ~isempty(con_f2f)
                ind_F = reduce(vcat,con_f2f)
                i = Vector{Int64}()
                for ii ∈  ind_F
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

function distmarch(F,V,indStart; d=missing, dd=missing, dist_tol=1e-3,con_V2V=missing,l=missing)

    # Get vertex-vertex connectivity
    if ismissing(con_V2V)
        con_V2V = con_vertex_vertex_f(F,V) 
    end

    # Compute "Laplacian umbrella" distances
    if ismissing(dd)
        dd = Dict{Vector{Int64},Float64}()  
        for i ∈ eachindex(V)
            for ii ∈ con_V2V[i]
                k = sort([i,ii])
                if !haskey(dd,k)
                    dd[sort(k)] = norm(V[i]-V[ii])
                end 
            end
        end
    end

    # Get/allocate distance vector
    if ismissing(d)
        d = fill(Inf,length(V))
    end

    if ismissing(l)
        l = fill(0,length(V))
    end

    # Set start distances to zero 
    d[indStart] .= 0.0
    l[indStart] .= 1:1:length(indStart)
    
    c = false
    ds = -1.0 # Set negative initially 

    boolCheck = fill(true,length(V))
    while true                          
        for i ∈ eachindex(V) # For each point            
            for ii ∈ con_V2V[i] # Check umbrella neighbourhood
                minVal,minInd = findmin([d[ii],dd[sort([i,ii])]+d[i]])            
                if minInd==2
                    d[ii] = minVal                          
                    l[ii] = l[i]
                end
            end            
        end
        if ~any(isinf.(d)) # Start checking once all are no longer Inf
            if c # If we were here before
                if abs(sum(d)-ds)<dist_tol                                        
                    break                    
                end
            end
            c = true # Flip to denote we've been here           
            ds = sum(d) # Now start computing the sum to check convergence
        end
    end
    return d,dd,l
end

# function distseedpoints(F,V,numPoints; ind=[1],dist_tol=1e-3)
    
#     con_V2V = con_vertex_vertex_f(F,V) 
#     d,dd,l = distmarch(F,V,ind; dist_tol=dist_tol,con_V2V=con_V2V)

#     if numPoints>1
#         @showprogress 1 "<distseedpoints>: Seeding points..." for q ∈ 2:1:numPoints            
#             push!(ind,findmax(d)[2])
#             d,dd,l = distmarch(F,V,ind; dist_tol=dist_tol, dd=dd,d=d,con_V2V=con_V2V,l=l)        
#         end
#     end
#     return ind,d,l
# end


"""
    ray_triangle_intersect(f::TriangleFace{Int64},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))

# Description 
This function can compute triangle-ray or triangle-line intersections. The 
function implements the "Möller-Trumbore triangle-ray intersection algorithm". 

# References 
[Möller, Tomas; Trumbore, Ben (1997). "Fast, Minimum Storage Ray-Triangle Intersection". Journal of Graphics Tools. 2: 21-28. doi: 10.1080/10867651.1997.10487468.](https://doi.org/10.1080/10867651.1997.10487468)
"""
function ray_triangle_intersect(f::TriangleFace{Int64},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))

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

    p = GeometryBasics.Point{3, Float64}(NaN,NaN,NaN)
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

function ray_triangle_intersect(F::Vector{TriangleFace{Int64}},V,ray_origin,ray_vector; rayType = :ray, triSide = 1, tolEps = eps(Float64))
    P = Vector{GeometryBasics.Point{3, Float64}}()
    indIntersect = Vector{Int64}()
    for qf ∈ eachindex(F)
        p = ray_triangle_intersect(F[qf],V,ray_origin,ray_vector; rayType = rayType, triSide = triSide, tolEps = tolEps)        
        if ~any(isnan.(p))
            push!(P,p)
            push!(indIntersect,qf)
        end
    end
    return P,indIntersect
end


"""
    mesh_curvature_polynomial(F,V)

# Description
Thisi function computes the per vertex curvature for the input mesh defined by 
the face `F` and the vertices `V`. A local polynomial is fitted (x+y+x^2+xy+y^2) to each points
"Laplacian umbrella" (point neighbourhood), and the curvature of this fitted
form is  

Implemented with the aid of [this helpful document](https://github.com/alecjacobson/geometry-processing-curvature/blob/master/README.md)

# References 
[F. Cazals and M. Pouget, "Estimating differential quantities using polynomial fitting of osculating jets", Computer Aided Geometric Design, vol. 22, no. 2, pp. 121-146, Feb. 2005, doi: 10.1016/j.cagd.2004.09.004](https://doi.org/10.1016/j.cagd.2004.09.004)
"""
function mesh_curvature_polynomial(F,V)
    
    E = meshedges(F)
    E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)  

    con_V2V = con_vertex_vertex(E_uni,V)


    N = vertexnormal(F,V)
    nz = Vec{3,Float64}(0.0,0.0,1.0)

    K1 = Vector{Float64}(undef,length(V))
    K2 = Vector{Float64}(undef,length(V))
    U1 = Vector{Vec3{Float64}}(undef,length(V))
    U2 = Vector{Vec3{Float64}}(undef,length(V))
    for q ∈ eachindex(V)
        n = N[q]
        Q = rotation_between(n,nz)
        ind = con_V2V[q]        
        vr = [Q*(v-V[q]) for v in V[ind]]

        T = Matrix{Float64}(undef,(length(ind),5))
        w = Matrix{Float64}(undef,(length(ind),1))
        for i = 1:1:length(ind)
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
        d,u = eigen(S) # Eigen decomposition to get first/second eigenvalue and vectors
        
        K1[q] = d[2]
        K2[q] = d[1]
        U1[q] = Q'*Vec3{Float64}(u[1,2],u[2,2],0.0)
        U2[q] = Q'*Vec3{Float64}(u[1,1],u[2,1],0.0)    
    end

    H = 0.5 * (K1.+K2)
    G = K1.*K2

    return K1,K2,U1,U2,H,G
end
    
function mesh_curvature_polynomial(M::GeometryBasics.Mesh) 
        return mesh_curvature_polynomial(faces(M),coordinates(V))
end

function seperate_vertices(F,V)
    Vn = Vector{eltype(V)}()
    Fn = deepcopy(F)
    c = 0 
    for q ∈ eachindex(F)
        f = F[q]
        m = length(f)
        Fn[q] = c .+ (1:m)
        c += m
        append!(Vn,V[f])
    end
    return Fn,Vn
end

function seperate_vertices(M::GeometryBasics.Mesh)
    F = faces(M)
    V = coordinates(M)
    Fn,Vn = seperate_vertices(F,V)    
    return GeometryBasics.Mesh(Vn,Fn)
end

function curve_length(V)
    return pushfirst!(cumsum(norm.(diff(V))),0.0) # Along curve distance from start
end

function evenly_sample(V, n; rtol = 1e-8, niter = 1)
    m = length(V)
    T = curve_length(V) # Initialise as along curve (multi-linear) distance
    T ./= last(T) # Normalise
    S = BSplineKit.interpolate(T, V, BSplineOrder(4), BSplineKit.Natural()) # Create interpolator
    for _ ∈ 1:niter
        dS = Derivative() * S  # spline derivative
        L = similar(T) # Initialise spline length vector
        L[1] = 0
        for i ∈ 2:m
            # Compute length of segment [i-1, i]
            segment_length, _ = quadgk(T[i-1], T[i]; rtol) do t
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
