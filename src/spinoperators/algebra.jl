function *(number::Real, ps::PauliString{S,N}) where {S,N}
    block = Block{S,N}([Float64(number)], [ps])
    return block
end

*(ps::PauliString{S,N}, number::Real) where {S,N} = *(number, ps)

function *(number::Real, block::Block{S,N}) where {S,N}
    new_block = Block{S,N}(Float64(number)*block.prefactors, block.pauli_strings)
    return new_block
end

*(block::Block, number::Real) = *(number, block)

function +(ps1::PauliString{S,N}, ps2::PauliString{S,N}) where {S,N}
    block = Block{S,N}([1.0, 1.0], [ps1, ps2])
    return block
end

function +(block::Block{S,N}, ps::PauliString{S,N}) where {S,N}
    new_block = Block{S,N}(vcat(block.prefactors, 1.0), vcat(block.pauli_strings, ps))
    return new_block
end

function +(ps::PauliString{S,N}, block::Block{S,N}) where {S,N}
    new_block = Block{S,N}(vcat(1.0, block.prefactors), vcat(ps, block.pauli_strings))
    return new_block
end

function +(block1::Block{S,N}, block2::Block{S,N}) where {S,N}
    new_block = Block{S,N}(
        vcat(block1.prefactors, block2.prefactors),
        vcat(block1.pauli_strings, block2.pauli_strings),
    )
    return new_block
end

function -(ps1::PauliString{S,N}, ps2::PauliString{S,N}) where {S,N}
    block = Block{S,N}([1.0, -1.0], [ps1, ps2])
    return block
end

function -(block::Block{S,N}, ps::PauliString{S,N}) where {S,N}
    new_block = Block{S,N}(vcat(block.prefactors, -1.0), vcat(block.pauli_strings, ps))
    return new_block
end

function -(ps::PauliString{S,N}, block::Block{S,N}) where {S,N}
    new_block = Block{S,N}(
        vcat(1.0, -1.0 * block.prefactors), vcat(ps, block.pauli_strings)
    )
    return new_block
end

function -(block1::Block{S,N}, block2::Block{S,N}) where {S,N}
    new_block = Block{S,N}(
        vcat(block1.prefactors, -1.0 * block2.prefactors),
        vcat(block1.pauli_strings, block2.pauli_strings),
    )
    return new_block
end

function ^(ps::PauliString{S,N,L}, pow::Int) where {S,N,L}
    return PauliString{S,N,L}(ps.sig, ps.locs, ps.power*pow)
end
