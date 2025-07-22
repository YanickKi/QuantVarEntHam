function subscriptnumber(i::Int)
    if i < 0
        c = [Char(0x208B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        push!(c, Char(0x2080+d))
    end
    return join(c)
end

function superscriptnumber(i::Int)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
   return join(c)
end

function float_to_int(x::Float64)
    isinteger(x) ? Int(x) : x
end

function print_pauli_string(io::IO, ps::PauliString)
    ps.power > 1 ?  print(io, "(") : nothing 
    for loc in ps.locs
        print(io, ps.sig * subscriptnumber(loc))
        last(ps.locs) == loc ? nothing : print(io, "⊗ ")
    end
   
    ps.power > 1 ? print(io, ")" * superscriptnumber(ps.power)) : nothing
end 

function print_block(io::IO, block::Block)

    if block.prefactors[1] < 0 
        print(io, "-")
    end 

    abs(block.prefactors[1]) == 1 ? nothing : print(io, float_to_int(abs(block.prefactors[1])), "*")

    print_pauli_string(io, block.pauli_strings[1])

    for index in 2:lastindex(block.pauli_strings) 
        block.prefactors[index] < 0 ? print(io, " - ") : print(io, " + ")
        abs(block.prefactors[index]) == 1 ? nothing : print(io, float_to_int(abs(block.prefactors[index])), "*")
        print_pauli_string(io,block.pauli_strings[index])
    end
end


function print_info(io::IO, ::AbstractBlock{S,N}) where {S,N}
    println(io ,"Spin " , S)
    println(io ,"Number of spins: ", N)
end 

function Base.show(io::IO, ::MIME"text/plain", ps::PauliString{S, N}) where {S,N}
    println(io, "Pauli string")
    print_info(io, ps)
    println(io)
    print_pauli_string(io, ps)
end

function Base.show(io::IO, ::MIME"text/plain", block::Block{S, N}) where {S,N}
    println(io, "Block")
    print_info(io, block)
    println(io)
    print_block(io, block)
end

function Base.show(io::IO, ::MIME"text/plain", pauli_strings::Vector{PauliString{S,N,L}}) where {S,N,L}
    println(io, "Number of Pauli strings: ", length(pauli_strings))
    println(io)
    print_info(io, pauli_strings[1])
    println(io)
    for index in eachindex(pauli_strings)
        println("Pauli string ", index, ": " ) 
        print(io, "\t \t")
        print_pauli_string(io, pauli_strings[index])
        if index != lastindex(pauli_strings)
            println(io)
            println(io)
        end 
    end 
end 

function Base.show(io::IO, ::MIME"text/plain", pauli_strings::Vector{PauliString{S,N}}) where {S,N}
    println(io, "Number of Pauli strings: ", length(pauli_strings))
    println(io)
    print_info(io, pauli_strings[1])
    println(io)
    for index in eachindex(pauli_strings)
        println("Pauli string ", index, ": " ) 
        print(io, "\t \t")
        print_pauli_string(io, pauli_strings[index])
        if index != lastindex(pauli_strings)
            println(io)
            println(io)
        end 
    end 
end 

function Base.show(io::IO, ::MIME"text/plain", blocks::Vector{Block{S,N}}) where {S,N}
    println(io, "Number of blocks: ", length(blocks))
    println(io)
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

Base.show(io::IO, block::AbstractBlock) = show(io, MIME"text/plain"(), block)

Base.show(io::IO, blocks::Vector{<:AbstractBlock}) = show(io, MIME"text/plain"(), blocks)
