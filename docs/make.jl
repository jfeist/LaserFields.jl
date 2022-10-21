using LaserFields
using Documenter

DocMeta.setdocmeta!(LaserFields, :DocTestSetup, :(using LaserFields); recursive=true)

makedocs(;
    modules=[LaserFields],
    authors="Johannes Feist <johannes.feist@gmail.com> and contributors",
    repo="https://github.com/jfeist/LaserFields.jl/blob/{commit}{path}#{line}",
    sitename="LaserFields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jfeist.github.io/LaserFields.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jfeist/LaserFields.jl",
    devbranch="main",
)
