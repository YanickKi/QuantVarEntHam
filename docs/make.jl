push!(LOAD_PATH,"../src/")

using Documenter, QuantVarEntHam, Yao

makedocs(sitename="QuantVarEntHam.jl",  modules = [QuantVarEntHam], 
        repo = "https://gitlab.dlr.de/ma-kind/QuantVarEntHam.git")