using GeneFinder
using Documenter

DocMeta.setdocmeta!(GeneFinder, :DocTestSetup, :(using GeneFinder); recursive=true)

makedocs(;
    modules=[GeneFinder],
    authors="Camilo García",
    repo="https://github.com/camilogarciabotero/GeneFinder.jl/blob/{commit}{path}#{line}",
    sitename="GeneFinder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://camilogarciabotero.github.io/GeneFinder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "A first algorithm" => "simplefinder.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/camilogarciabotero/GeneFinder.jl",
    devbranch="main",
)
