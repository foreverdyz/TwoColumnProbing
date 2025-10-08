#coupling_matrix.jl

using SparseArrays

"""
    Now we need to build coupling matrix: c_mt. c_mt[i,j] implies the number of constraints containing the i-th and j-th binary
variables simultaneously.
    We will check all binary set from B.

    coupling_matrix()
This function builds a sparse matrix for c_mt (aforementioned).
"""
function coupling_matrix(B::Vector{Vector{Int64}}, binary_number::Int64, size_limit::Int64, work_limit::Int64)
    let
        local a = Int64[]
        local b = Int64[]
        local s = 0
        for set in B
            local st = 2
            local en = length(set)
            local id = true
            #ignore too long constraints because the memory cost
            (en > size_limit) && (id = false)
            if id
                for i in st : en
                    #index-1 to avoid diag term
                    for k in 1:(i-1)
                        if set[i] > set[k]
                            push!(a, set[i])
                            push!(b, set[k])
                        else
                            push!(a, set[k])
                            push!(b, set[i])
                        end
                        s += 1
                    end
                end
            end
            #that is the totally workload
            (s > work_limit) && (break)
        end
        push!(a, binary_number)
        push!(b, binary_number)
        local c = ones(Int64, s)
        push!(c, 0)
        return dropzeros!(sparse(a, b, c))
    end
end