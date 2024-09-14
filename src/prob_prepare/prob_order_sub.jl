#prob_order_sub.jl

using SparseArrays
using BangBang

"""
    Now we need to sort all pair of binary variables with some polices:
    first, the most important factor is the number of constraints that two variables appearing simultaneously;
    second, another factor is the number of variables conflicts the two variables;
    third, we do not want the two variables have conflicted; thus if they have 1 conflict, the score is divided by 2, and if there
is 2 conflicts, the two variables will be fixed by each other (and we ignore them);
    fourth, we should also consider about the implementation graph. However, since we cannot use SCIP to presolve at first, we give
up this factor in this work.
"""

function prob_order_sub(
        con_matrix::AbstractSparseArray, clq_table::Vector{Vector{Int64}}, clq_length::Vector{Int64}, c_mt::AbstractSparseArray, 
        binary_length::Int64, assign_list_sub::Vector{Int64}, threadnum::Int64, bin_to_org::Dict{Int64, Int64},
        cand_number::Int64, max_probe_number::Int64
    )
    let
        local I, J, K = findnz(c_mt[assign_list_sub, assign_list_sub])
        local n = ceil(Int64, cand_number/log(2, threadnum+1))
        if length(I) > n
            order = _select_candidate(K, n)
            #order = partialsortperm(K, rev = true, 1:n)
            I, J, K = I[order], J[order], K[order]
        end
        I, J = assign_list_sub[I], assign_list_sub[J]
        
        var_con_len = zeros(binary_length)
        for i in 1:binary_length
            x = bin_to_org[i]
            var_con_len[i] = nnz(con_matrix[:, x])
        end
        for i in 1:length(I)
            local conflict = 1
            local x = I[i]
            local y = J[i]
            (_have_conflict(clq_table[x], clq_table[y]) > 0) && (conflict = 10)
            K[i] = ceil(Int64, (10*K[i] + var_con_len[x] + var_con_len[y] + 3*clq_length[x] + 3*clq_length[y])/conflict)
        end
        
        local m = min(length(I), ceil(Int64, max_probe_number/log(2, threadnum+1)))
        local counter_order = partialsortperm(K, rev = true, 1:m)
        return I[counter_order], J[counter_order]
    end
end


function _have_conflict(tb_x::Vector{Int64}, tb_y::Vector{Int64})
    st = 1
    en = length(tb_y)
    for x in tb_x
        for i in st:en
            if x == tb_y[i]
                return 1
            else
                if x < tb_y[i]
                    st = i
                    break
                end
            end
        end
    end
    return 0
end

function _select_candidate(K::Vector{Int64}, n::Int64)
    local upper = findmax(K)[1]
    local step = (upper - 1)/2
    local lower = 1
    local pivot = ceil(Int64, 1 + step)
    local I = [Int64(i) for i in 1:length(K)]
    local p = Int64[]
    while pivot <= upper
        push!(p, pivot)
        local U = Int64[]
        local L = Int64[]
        local lower_id = false
        for i in I
            (K[i] >= pivot) ? (push!(U, i)) : (push!(L, i))
        end
        local m = length(U)
        step = max(step/2, 1)
        if m > n
            lower = max(pivot, lower)
            pivot = ceil(Int64, pivot + step)
            I = copy(U)
        else
            if 1.5*m < n
                pivot = ceil(Int64, lower + step)
                lower_id = true
            else
                return U
            end
        end
        if pivot in p
           return U
        end
        (lower_id) && (I = L)
    end
    return I
end