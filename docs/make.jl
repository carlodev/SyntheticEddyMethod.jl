push!(LOAD_PATH,"../src/")
using SyntheticEddyMethod
using Documenter
makedocs(
         sitename = "SyntheticEddyMethod.jl",
         modules  = [SyntheticEddyMethod],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/carlodev/SyntheticEddyMethod.jl",)