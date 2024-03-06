using Documenter, Comodo

makedocs(
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  collapselevel = 2,
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
