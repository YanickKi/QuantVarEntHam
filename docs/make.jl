using Documenter, QuantVarEntHam

makedocs(sitename="QuantVarEntHam.jl",  
        modules = [QuantVarEntHam], 
        pages = [
                "Home" => "index.md"
                "Tutorials" => "tutorials.md"
                "Spinoperators" => "spinoperators.md"
                "Models" => "models.md"
                "Ansätze" => "ansätze.md"
                "Cost" => "cost.md"
                "Integration" => "integration.md"
                "Optimizer" => "optimizer.md"
        ],
        repo = Documenter.Remotes.GitLab("gitlab.dlr.de", "ma-kind", "QuantVarEntHam"))