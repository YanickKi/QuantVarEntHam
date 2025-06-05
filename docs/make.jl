push!(LOAD_PATH,"../src/")

using Documenter, QuantVarEntHam

makedocs(sitename="QuantVarEntHam.jl",  modules = [QuantVarEntHam], 
        pages = [
                "Home" => "index.md"
                "Functionality" => "functionality.md"
                "Models" => "models.md"
                "Initialize" => "initialize.md"
                "Ansätze" => "ansätze.md"
                "Optimizer" => "optimizer.md"
        ],
        repo = "https://gitlab.dlr.de/ma-kind/QuantVarEntHam.git")