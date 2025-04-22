using Yao

@with_kw mutable struct Settings_toric{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T, S}
    Nx::Int
    Ny::Int
    A::Vector{Int}
    N_A::Int
    J::Float64
    T_max::Float64
    ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} 
    observables::Vector{T}
    meas0::Vector{Float64}  = [expect(observables[i], ρ_A) for i in eachindex(observables)]
    mtrxObs::Vector{S}
end

function toric(Nx::Int, Ny::Int, A::Vector{Int}, J::Real, T_max::Real; ρ_A::DensityMatrix{2}=get_rhoA(H_toric(Nx, Ny, J),  A, Nx ÷ 2 * Ny),
    observables::Vector{<:AbstractBlock}=[put(length(A), qbit => X) for qbit in eachindex(A)])

    mtrxObs = mat.(observables)

    return Settings_toric{eltype(observables), eltype(mtrxObs)}(
        Nx = Nx, Ny = Ny, A = A, N_A = length(A), J = J, T_max = T_max,
        ρ_A = ρ_A, observables = observables,
        mtrxObs = mtrxObs,
    ) 
end 


function map_to_chain(obj, Nx)
    qbitnumbers = Int64[]
    _Nx = Nx ÷ 2

    for qbit in 1:4
        x =  obj[qbit][1]
        y =  obj[qbit][2]
        
        row = 0
        
        if x % 2 == 0
            row = 0
        else 
            row = 1 
        end 

        _x = (x+row) ÷ 2 

        push!(qbitnumbers, _Nx * (y-1) + _x)
    end 
    return qbitnumbers
end 

function make_all_objects(Nx,Ny)
    stars = Vector{Tuple{Int64, Int64}}[]
    plaquets = Vector{Tuple{Int64, Int64}}[]

    for y in 1:2:Ny-1
        for x in 2:2:Nx
        push!(stars, [(x,y), 
                    (x == 2 ? Nx : x-2, y), 
                    (x-1, y+1), 
                    (x-1, y == 1 ? Ny : y-1)])
        end 
    end

    for y in 2:2:Ny
        for x in 1:2:Nx-1
        push!(plaquets, [(x,y), 
                    (x == Nx-1 ? 1 : x+2, y), 
                    (x+1, y == Ny ? 1 : y+1), 
                    (x+1, y-1)])
        end 
    end


    stars_qubits = Vector{Int64}[]
    plaquets_qubits = Vector{Int64}[]
    for i in eachindex(stars)
        push!(stars_qubits, map_to_chain(stars[i], Nx))
        push!(plaquets_qubits, map_to_chain(plaquets[i], Nx))
    end

    return stars_qubits, plaquets_qubits

end 


function H_toric(Nx, Ny, J)
 
    stars_qubits, plaquets_qubits = make_all_objects(Nx,Ny)
    
    N_obj = length(stars_qubits)
    N_qubits = (Nx ÷ 2) * Ny

    H_s = map(1:N_obj) do i 
        repeat(N_qubits, X, stars_qubits[i])
    end  |> sum 
    
    H_p = map(1:N_obj) do i 
        repeat(N_qubits, Z, plaquets_qubits[i])
    end  |> sum 

    println(H_s+H_p)

    return -J * (H_s + H_p)

end 