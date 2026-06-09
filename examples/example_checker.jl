module ExampleChecker
using Comodo

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
    example_dir = joinpath(comododir(),"examples")
    dircontent = readdir(example_dir)
    juliafiles = filter(x -> occursin(r"demo.*\.jl", x), dircontent)
    
    problems = Vector{String}(undef, 0)    
    for file in juliafiles        
        println("Checking file: $file")
        t = Task(() -> load_and_run!(joinpath(example_dir, file), problems))
        schedule(t)        
        while !istaskdone(t)  # Wait for task completion                  
            sleep(1)
        end
    end

    if !isempty(problems)
        @info "Here is the list of problematic files:"
        display(problems)
    else
        @info "All demos completed successfully."
    end
end

function __init__()
    @info "Module initialized, running main"
    main()
end

end # end of module
