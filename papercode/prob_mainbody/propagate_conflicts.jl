#propagate_conflicts.jl

using SparseArrays

"""
    propagate_conflicts()
Fixing variables based on clique table. For example, fixing x = 1, we may fix y = 0 since x + y <= 1.
"""
function propagate_conflicts(
        x::Int64, is_p::Bool,
        clq_table::Vector{Vector{Int64}}, clq_p_n::Vector{Vector{Int64}},
        S::Vector{Vector{Int64}}, binary_number::Int64
    )
    conflict, has_prop = spzeros(binary_number), spzeros(binary_number)
    wait_p, wait_n = Int64[], Int64[]
    if is_p
        push!(wait_p, x)
        conflict[x] = 1
    else
        push!(wait_n, x)
        conflict[x] = -1
    end
    while length(wait_p) + length(wait_n) > 0
        if length(wait_p) > 0
            y = popfirst!(wait_p)
            con_list = clq_table[y]
            for (i, j) in enumerate(con_list)
                if clq_p_n[y][i] > 0
                    for k in S[j]
                        if k <= binary_number && k != y
                            #k will be fixed to 0
                            (conflict[k] > 0) ? (return conflict, true) : (conflict[k] = -1)
                            (has_prop[k] < 1) && (push!(wait_n, k))
                        elseif k > binary_number
                            #k - binary_number will be fixed to 1
                            (conflict[k - binary_number] < 0) ? (return conflict, true) : (conflict[k - binary_number] = 1)
                            (has_prop[k - binary_number] < 1) && (push!(wait_p, k - binary_number))
                        end
                    end
                end
            end
            has_prop[y] = 1
        end
        if length(wait_n) > 0
            y = popfirst!(wait_n)
            con_list = clq_table[y]
            for (i, j) in enumerate(con_list)
                if clq_p_n[y][i] < 0
                    for k in S[j]
                        if k <= binary_number
                            #k will be fixed to 0
                            (conflict[k] > 0) && (return conflict, true) : (conflict[k] = -1)
                            (has_prop[k] < 1) && (push!(wait_n, k))
                        elseif k > binary_number && k != y + binary_number
                            #k - binary_number will be fixed to 1
                            (conflict[k - binary_number] < 0) ? (return conflict, true) : (conflict[k - binary_number] = 1)
                            (has_prop[k - binary_number] < 1) && (push!(wait_p, k - binary_number))
                        end
                    end
                end
            end
            has_prop[y] = 1
        end
    end
    return conflict, false
end