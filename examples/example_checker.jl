module ExampleChecker

function load_and_run!(filename::String, problems::Vector{String})::Bool
    cmdline = "using Pkg; Pkg.activate(\"./..\"); include(\"$filename\")"
    command = `julia -e $cmdline`
    try
        result =  run(command)
        if result.exitcode != 0
            push!(problems, filename)
            return false
        end
        return true
    catch e
        push!(problems, filename)
        return false
    end
end

function main()
    dircontent = readdir(".")
    juliafiles = filter(x -> occursin(r"demo.*\.jl", x), dircontent)

    problems = Vector{String}(undef, 0)
    for file in juliafiles
        println("Checking file: $file")
        load_and_run!(file, problems)
    end

    if !isempty(problems)
        @info "Here is the list of problematic files:"
        display(problems)
    else
        @info "All demos completed succesfully."
    end
end

function __init__()
    @info "Module initialized, running main"
    main()
end

end # end of module