using Documenter, QuantVarEntHam

makedocs(sitename="QuantVarEntHam.jl",  
        modules = [QuantVarEntHam], 
        pages = [
                "Home" => "index.md"
                "Tutorials" => "tutorials.md"
                "Spin Operators" => "spinoperators.md"
                "Models" => "models.md"
                "Ansätze" => "ansätze.md"
                "Cost" => "cost.md"
                "Optimizer" => "optimizer.md"
                "Universal Ratios" => "universal_ratios.md"
        ],
        repo = Documenter.Remotes.GitLab("gitlab.dlr.de", "ma-kind", "QuantVarEntHam"))