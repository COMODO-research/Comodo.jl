# Usage:
# > julia coverage.jl 
# in terminal

import Coverage 
import Pkg 

Pkg.activate("../")
Pkg.resolve()

Pkg.test(coverage = true)

tested, total = Coverage.process_folder("../") |> Coverage.get_summary

ratio = tested / total 

Coverage.clean_folder("../")


@info "Coverate: $ratio"