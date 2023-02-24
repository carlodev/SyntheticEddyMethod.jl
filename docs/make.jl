using SyntheticEddyMethod, Documenter
makedocs(sitename="SyntheticEddyMethod")

makedocs(
    strict=true,
    sitename = "SyntheticEddyMethod.jl",
    modules = [PartitionedArrays],
    pages = [
        "Introduction" => "index.md",
        "usage.md",
        "exploring.md",
    ],
)
