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

deploydocs(
    repo = "github.com/carlodev/SyntheticEddyMethod.jl.git",
    push_preview = true,
)
