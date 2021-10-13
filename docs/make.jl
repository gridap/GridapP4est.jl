using GridapP4est
using Documenter

DocMeta.setdocmeta!(GridapP4est, :DocTestSetup, :(using GridapP4est); recursive=true)

makedocs(;
    modules=[GridapP4est],
    authors="Alberto F. Martin <alberto.martin@monash.edu>",
    repo="https://github.com/gridap/GridapP4est.jl/blob/{commit}{path}#{line}",
    sitename="GridapP4est.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gridap.github.io/GridapP4est.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gridap/GridapP4est.jl",
    devbranch="main",
)
