module ExampleChecker

const NUMTASKS = 10

function load_and_run!(filename::String, problems::Vector{String})::Bool
    #cmdline = "using Pkg; Pkg.activate(\"./..\"); include(\"$filename\")"
		cmdline = "using Pkg; Pkg.activate(\".\"); Pkg.instantiate(); include(\"$filename\")"
    command = `julia --project=. -e $cmdline`
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
    tasks = Task[]
    for file in juliafiles
        if length(tasks) >= NUMTASKS
            wait.(tasks)
            empty!(tasks)
        end
        println("Checking file: $file")
        t = Task(() -> load_and_run!(file, problems))
        push!(tasks, t)
        schedule(t)
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
