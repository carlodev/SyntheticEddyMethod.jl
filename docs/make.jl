using  Documenter, SyntheticEddyMethod
makedocs(
    sitename = "SyntheticEddyMethod.jl",
    modules = [SyntheticEddyMethod],
    pages = [
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "Advanced Usage" => "advanced_usage.md",
        "Exploring" => "exploring.md",
    ],
)
