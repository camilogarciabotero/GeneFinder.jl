using GeneFinder
using Documenter

DocMeta.setdocmeta!(GeneFinder, :DocTestSetup, :(using GeneFinder); recursive = true)

makedocs(;
    modules = [GeneFinder],
    authors = "Camilo García",
    repo = "https://github.com/camilogarciabotero/GeneFinder.jl/blob/{commit}{path}#{line}",
    sitename = "GeneFinder.jl",
    format = Documenter.HTML(;
        mathengine=MathJax3(),
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://camilogarciabotero.github.io/GeneFinder.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Finding ORFs" => "simplefinder.md",
        "Towards Markov Chains" => "markovchains.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/camilogarciabotero/GeneFinder.jl", devbranch = "main")
