# Run with: julia --project=./docs/make.jl
using Documenter, DocStringExtensions, GCATBase
makedocs(modules = [GCATBase], sitename="GCATBase.jl"; remotes = nothing)