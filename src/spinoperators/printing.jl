function print_pauli_string(io::IO, ps::PauliString)
    for _ in 1:ps.power
        for loc in ps.locs
            print(io, ps.sig, loc)
            last(ps.locs) == loc ? nothing : print(io, *)
        end
    end 
end 

function print_block(io::IO, block::Block)
    for index in eachindex(block.pauli_strings) 
        block.prefactors[index] == 1 ? nothing : print(io, block.prefactors[index], " * ")
        print_pauli_string(io,block.pauli_strings[index])
        index == length(block.pauli_strings) ? nothing : print(io, " + ")
    end
end

function print_info(io::IO, ::AbstractBlock{S,N}) where {S,N}
    println(io ,"Spin " , S)
    println(io ,"Number of spins: ", N)
end 

function Base.show(io::IO, ps::PauliString{S, N}) where {S,N}
    print_info(io, ps)
    println(io)
    print_pauli_string(io, ps)
    println(io)
end

function Base.show(io::IO, block::Block{S, N}) where {S,N}
    print_info(io, block)
    println(io, )
    print_block(io, block)
end

function Base.show(io::IO, ::MIME"text/plain", blocks::Vector{Block{S,N}}) where {S,N}
    println(io, "Number of blocks: ", length(blocks))
    print_info(io, blocks[1])
    println(io)
    for index in eachindex(blocks)
        println("Block ", index, ": " ) 
        print(io, "\t")
        print_block(io, blocks[index])
        if index != lastindex(blocks)
            println(io)
            println(io)
        end 
    end 
end 