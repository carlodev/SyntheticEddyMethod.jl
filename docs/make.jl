using  Documenter, SyntheticEddyMethod
makedocs(
    sitename = "SyntheticEddyMethod.jl",
    modules = [SyntheticEddyMethod],
    pages = [
        "Introduction" => "index.md",
        "Usage" => "usage.md",
        "Exploring" => "exploring.md",
    ],
)
