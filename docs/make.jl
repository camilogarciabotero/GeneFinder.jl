using GeneFinder
using Documenter

DocMeta.setdocmeta!(GeneFinder, :DocTestSetup, :(using GeneFinder); recursive = true)

makedocs(;
    modules = [GeneFinder],
    authors = "Camilo García",
    repo = "https://github.com/camilogarciabotero/GeneFinder.jl/",
    sitename = "GeneFinder.jl",
    format = Documenter.HTML(;
        mathengine=MathJax3(),
        # prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://camilogarciabotero.github.io/GeneFinder.jl",
        repolink = "https://github.com/camilogarciabotero/GeneFinder.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Finding ORFs" => "naivefinder.md",
        # "A Simple coding rule" => "simplecodingrule.md",
        "Wrtiting ORFs in files" => "iodocs.md",
        "Roadmap" => "roadmap.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "https://github.com/camilogarciabotero/GeneFinder.jl", devbranch = "main")