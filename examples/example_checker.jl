module ExampleCheckter

using Pkg
Pkg.activate("../.")

using Test

struct Problem
    filename::String
    error::Exception
end

function Base.show(io::IO, p::Problem)
    println(io, "Problem in file: $(p.filename)")
    println(io, "Error: $(p.error)")
end

function load_and_run!(filename::String, problems::Vector{Problem})::Bool
    try
        @test begin
            include(filename)
            return true
        end
    catch e
        p = Problem(filename, e)
        push!(problems, p)
        println("Error in file: $filename")
        println("Please check the error message to find out the required function and package")
        println(e)
        return false
    end
end

function main()
    dircontent = readdir(".")
    juliafiles = filter(x -> occursin(r"demo.*\.jl", x), dircontent)

    problems = Vector{Problem}(undef, 0)

    for file in juliafiles
        println("Checking file: $file")
        load_and_run!(file, problems)
    end

    display(problems)
end

function __init__()
    @info "Module initialized, running main"
    main()
end

end # end of module