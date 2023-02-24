using SyntheticEddyMethod, Documenter
makedocs(sitename="SyntheticEddyMethod")

makedocs(
    strict=true,
    sitename = "SyntheticEddyMethod.jl",
    pages = [
        "Introduction" => "index.md",
        "usage.md",
        "exploring.md",
    ],
)
