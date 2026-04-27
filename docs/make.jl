# using Documenter, Comodo

# makedocs(
#          format = Documenter.HTML(
#                                   prettyurls = get(ENV, "CI", nothing) == "true",
#                                   collapselevel = 2,
#                 		  size_threshold_warn=1500 * 2^10, 
#         			  size_threshold=2000 * 2^10, 
#                                   # assets = ["assets/favicon.ico", "assets/extra_styles.css"],
#                                  ),
#          sitename="A Julia package for computational (bio)mechanics and computational design",
#          authors = "Kevin-Mattheus-Moerman <kevin.moerman@gmail.com>",
#          pages = [
#                   "Functions" => "functions.md",
#                  ]
#         )


# deploydocs(
#            repo = "https://github.com/COMODO-research/Comodo.jl",
#           )


# using Comodo
# using Documenter
# using DocumenterVitepress

# DocMeta.setdocmeta!(Comodo, :DocTestSetup, :(using Comodo); recursive=true)

# makedocs(;
#     modules=[Comodo],
#     authors= "Kevin-Mattheus-Moerman <kevin.moerman@gmail.com>",
#     sitename="A Julia package for computational (bio)mechanics and computational design",
#     format=DocumenterVitepress.MarkdownVitepress(
#         repo = "github.com/Aminofa70/Comodo.jl",
#         devbranch="main",
#         devurl="dev",
#     ),
#     pages = [
#     "Home" => "index.md",

#     "Tutorials" => [
#         "Getting started" => "tutorials/getting-started.md",

#         "Demos" => [
#             "Demo 0001" => "tutorials/demo_0001.md",
#            # "Demo 0002" => "tutorials/demo_0002.md",
#         ],
#     ],

#     "API Reference" => "api.md",
# ]
# )

# DocumenterVitepress.deploydocs(;
#     repo = "github.com/Aminofa70/Comodo.jl.git",
#     target = joinpath(@__DIR__, "build"),
#     branch = "gh-pages",
#     devbranch = "main",
#     push_preview = true,
# )


using Comodo
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(Comodo, :DocTestSetup, :(using Comodo); recursive=true)

makedocs(;
    modules = Module[],
    authors = "Kevin-Mattheus-Moerman <kevin.moerman@gmail.com>",
    sitename = "Comodo.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/COMODO-research/Comodo.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting started" => "tutorials/getting-started.md",
            "Demos" => [
                "demo_quadplate" => "tutorials/demo_quadplate.md",
            ],
        ],
        "API Reference" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/COMODO-research/Comodo.jl.git",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)