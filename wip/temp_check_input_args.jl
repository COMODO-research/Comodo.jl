
function checkargset(args...)
    # The argument is a Tuple, e.g. Tuple{Int64, Irrational{:π}, String}
    if isnothing(args)
        args = ""
    else
        for q ∈ eachindex(args)
            a = args[q]
            # if isnothing(a)
            #     args[q] = "" # Empty string
            # else isa(args[q],Irrational) # Convert irrationals to numerical representation to avoid name being printed π
            #     args[q] = Float64.(args[q])
            # end
            display(typeof(a))
        end        
    end
    return args
end

a = 1
b = pi 
c = "yes"
d = zeros(3,3)
argSet=checkargset(a,b,c,d,nothing)