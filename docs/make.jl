using LaserFields
using Documenter

DocMeta.setdocmeta!(LaserFields, :DocTestSetup, :(using LaserFields); recursive=true)

makedocs(;
    modules=[LaserFields],
    authors="Johannes Feist <johannes.feist@gmail.com> and contributors",
    sitename="LaserFields.jl",
    format=Documenter.HTML(;
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
