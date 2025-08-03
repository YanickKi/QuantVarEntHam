# COV_EXCL_START
"""
    getrho_A(model::AbstractModel)

Return the reduced density matrix of the `model`.
"""
getrho_A(model::AbstractModel) = model.ρ_A

"""
    getblocks(ansatz::AbstractAnsatz)

Return the blocks of the given `ansatz`.
"""
getblocks(ansatz::AbstractAnsatz) = ansatz.blocks

"""
    getr_max(ansatz::AbstractAnsatz)

Return the maximum range of interaction of the given `ansatz`.
"""
getr_max(ansatz::AbstractAnsatz) = ansatz.r_max
# COV_EXCL_STOP
