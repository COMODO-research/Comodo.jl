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

function _interp_biharmonic_spline(x::Union{Vector{T}, AbstractRange{T}},y::Union{Vector{T}, AbstractRange{T}},xi::Union{Vector{T}, AbstractRange{T}}; extrapolate_method=:linear,pad_data=:linear) where T<:Real

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
        throw(ArgumentError("Invalid pad_data method provided, valid options are :linear, :constant, and :none"))
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
        error("InvalidParameter: Invalid extrapolate_method method provided, valid options are :linear, :constant, and :biharmonic")
    end

    return yi
end