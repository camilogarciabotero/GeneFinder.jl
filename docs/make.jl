using GeneFinder
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(GeneFinder, :DocTestSetup, :(using GeneFinder); recursive = true)

# deploydocs(; repo = "https://github.com/camilogarciabotero/GeneFinder.jl", devbranch = "main")

fmt = DocumenterVitepress.MarkdownVitepress(
    repo = "https://github.com/camilogarciabotero/GeneFinder.jl",
    devbranch = "main",
)

pgs = [
    "Home" => "index.md",
    "Get Started" => "getstarted.md",
    "Usage" => [
        "The ORF type" => "orftype.md",
        "Scoring ORFs" => "features.md",
        "A Simple Coding Rule" => "simplecodingrule.md",
        "Ribosome Binding Sites" => "rbs.md",
        "Writing ORFs In Files" => "iodocs.md",
    ],
    "API" => "api.md",
]

makedocs(;
    modules = [GeneFinder],
    authors = "Camilo García-Botero",
    repo = Remotes.GitHub("camilogarciabotero", "GeneFinder.jl"),
    sitename = "GeneFinder.jl",
    format = fmt,
    pages = pgs,
    warnonly = true,
)

DocumenterVitepress.deploydocs(; 
    repo = "https://github.com/camilogarciabotero/GeneFinder.jl",
    devbranch = "main",
    target = "build", # this is where Vitepress stores its output
    # branch = "gh-pages",
    push_preview = true
)