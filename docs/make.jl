push!(LOAD_PATH,"../src/")

using Documenter, QuantVarEntHam, Yao

makedocs(sitename="QuantVarEntHam.jl",  modules = [QuantVarEntHam], 
        pages = [
                "Home" => "index.md"
                "Functionality" => "functionality.md"
                "Settings" => "settings.md"
                "Initialize" => "initialize.md"
                "Ansätze" => "ansätze.md"
                "Hamiltonians" => "hamiltonians.md"
                "Optimizer" => "optimizer.md"
        ],
        repo = "https://gitlab.dlr.de/ma-kind/QuantVarEntHam.git")