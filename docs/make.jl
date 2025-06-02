using Documenter
using SymbolicFact

makedocs(
    sitename = "SymbolicFact",
    format = Documenter.HTML(),
    modules = [SymbolicFact]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "https://github.com/8Rudy8/SymbolicFact.jl"
)=#
