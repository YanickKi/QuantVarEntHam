"""
    getrho_A(model::AbstractModel)

Return the reduced density matrix of the `model`.
"""
getrho_A(model::AbstractModel) = model.ρ_A
