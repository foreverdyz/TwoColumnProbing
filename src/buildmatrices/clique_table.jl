#build_cg.jl

using SparseArrays

"""
    Now we need to build a clique table to detect conflicts between binary variables.
    How to use clique table (clq_table): for i-th term, it is for i-th binary variable (bin_to_org[i] in the original variable list)
clq_table[i] includes a list, e.g. [x1, x2, x3, ...], then you can find cliques containing i by S[x1], S[x2],... or con_matrix[J[x1]],...
"""
function clique_table(S::Vector{Vector{Int64}}, binary_number::Int64)
    #clq is a lenght(S) * binary_number sparse matrix, which cachce all cliques
    clq = spzeros(length(S), binary_number)
    for j in 1:length(S)
        local s = S[j]
        for i in s
            if i <= binary_number
                clq[j, i] = 1
            else
                clq[j, i - binary_number] = -1
            end
        end
    end
    #then cache clq with a clique table
    clq_table = Vector{Int64}[]
    clq_p_n = Vector{Int64}[]
    clq_length = Int64[]
    for i = 1:binary_number
        local con = clq[:, i]
        local x, y = findnz(con)
        local len = 0
        push!(clq_table, x)
        push!(clq_p_n, y)
        for j in x
            len += length(S[j])
        end
        push!(clq_length, len)
    end
    return clq_table, clq_p_n, clq_length
end