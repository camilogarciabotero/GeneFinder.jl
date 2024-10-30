using GeneFinder
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(GeneFinder, :DocTestSetup, :(using GeneFinder); recursive = true)

# makedocs(;
#     modules = [GeneFinder],
#     authors = "Camilo García-Botero",
#     repo = "https://github.com/camilogarciabotero/GeneFinder.jl/",
#     sitename = "GeneFinder.jl",
#     format = Documenter.HTML(;
#         mathengine=MathJax3(),
#         # prettyurls = get(ENV, "CI", "false") == "true",
#         canonical = "https://camilogarciabotero.github.io/GeneFinder.jl",
#         repolink = "https://github.com/camilogarciabotero/GeneFinder.jl",
#         edit_link = "main",
#     ),
#     pages = [
#         "Home" => "index.md",
#         "Finding ORFs" => "naivefinder.md",
#         "Scoring ORFs" => "features.md",
#         "A Simple Coding Rule" => "simplecodingrule.md",
#         "Wrtiting ORFs In Files" => "iodocs.md",
#         "API" => "api.md",
#         # "Roadmap" => "roadmap.md",
#     ],
#     warnonly = true,
# )

deploydocs(; repo = "https://github.com/camilogarciabotero/GeneFinder.jl", devbranch = "main")


fmt = DocumenterVitepress.MarkdownVitepress(
    repo = "https://github.com/camilogarciabotero/GeneFinder.jl",
    devbranch = "main",
)

pgs = [
    "Home" => "index.md",
    "Finding ORFs" => "naivefinder.md",
    "Scoring ORFs" => "features.md",
    "A Simple Coding Rule" => "simplecodingrule.md",
    "Wrtiting ORFs In Files" => "iodocs.md",
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

deploydocs(; 
    repo = "https://github.com/camilogarciabotero/GeneFinder.jl",
    devbranch = "main",
    target = "build", # this is where Vitepress stores its output
    branch = "gh-pages",
    push_preview = true
)