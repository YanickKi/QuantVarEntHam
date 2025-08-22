using Documenter, QuantVarEntHam
using DocumenterCitations
using Pkg


PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaDocs/DocumenterCitations.jl"

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)


makedocs(;
    format=Documenter.HTML(
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/DocumenterCitations.jl",
        assets=String["assets/citations.css"],
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)."
    ),
    sitename="QuantVarEntHam.jl",
    warnonly=:doctest,
    modules=[QuantVarEntHam],
    pages=[
        "Home" => "index.md"
        "Tutorials" => "tutorials.md"
        "Spin Operators" => "spinoperators.md"
        "Models" => "models.md"
        "Ansätze" => "ansätze.md"
        "Cost" => "cost.md"
        "Optimizer" => "optimizer.md"
        "Universal Ratios" => "universal_ratios.md"
    ],
    plugins=[bib],
    repo=Documenter.Remotes.GitLab("gitlab.com", "quantum-computing-software", "QuantVarEntHam"),
)
