using Documenter, DocumenterCitations, SyntheticEddyMethod


bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)


makedocs(bib;
    sitename = "SyntheticEddyMethod.jl",
    modules = [SyntheticEddyMethod],
    pages = [
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "Advanced Usage" => "advanced_usage.md",
        "Exploring" => "exploring.md",
        "API information" => "api_info.md",
        "References" => "references.md",
    ],
)

deploydocs(
    repo = "github.com/carlodev/SyntheticEddyMethod.jl",
    push_preview = true,
)
