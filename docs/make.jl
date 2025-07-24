using Documenter, Comodo

makedocs(
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  collapselevel = 2,
                		  size_threshold_warn=1500 * 2^10, 
        			  size_threshold=2000 * 2^10, 
                                  # assets = ["assets/favicon.ico", "assets/extra_styles.css"],
                                 ),
         sitename="A Julia package for computational (bio)mechanics and computational design",
         authors = "Kevin-Mattheus-Moerman <kevin.moerman@gmail.com>",
         pages = [
                  "Functions" => "functions.md",
                 ]
        )


deploydocs(
           repo = "https://github.com/COMODO-research/Comodo.jl",
          )
