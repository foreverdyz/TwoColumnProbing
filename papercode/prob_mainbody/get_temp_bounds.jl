#get_temp_bounds.jl

using SparseArrays

"""
    get_temp_bounds()
Conduct implication graph propagation based on fixed variables from the conflict graph (clique table here) propagation.
Note that:
    aggre[i]: [[j, i_p_j, i_n_j], ...]
    imp_p[i]: [[j, 1, j_ub] or [j, 0, j_lb], ...]
"""
function get_temp_bounds(
        conflict::AbstractVector, binary_number::Int64,
        A::AbstractSparseArray, var_ub::AbstractVector, var_lb::AbstractVector,
        aggre_list::AbstractVector, aggre::AbstractVector,
        imp_p_list::AbstractVector, imp_p::AbstractVector, imp_n_list::AbstractVector, imp_n::AbstractVector,
        two_imp_list::AbstractSparseArray, two_imp::AbstractVector, bin_to_org::Dict{Int64, Int64}
    )
    temp_ub, temp_lb = copy(var_ub), copy(var_lb)
    #construct temporary bounds
    I = Int64[]
    J = Int64[]
    for i in findnz(conflict)[1]
        if conflict[i] > 0
            push!(J, i)
            #change conflict bounds
            if temp_lb[bin_to_org[i]] < 1
                temp_lb[bin_to_org[i]] = 1
                push!(I, bin_to_org[i])
                (temp_ub[bin_to_org[i]] < 1) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
            end
            #for aggre
            if aggre_list[i] > 0
                record = aggre[Int64(aggre_list[i])]
                for rd in record
                    j = Int64(rd[1])
                    b = rd[2]
                    if temp_lb[j] > b || temp_ub[j] < b
                        return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]) #the last true implies we do not run this scenario
                    else
                        if temp_ub[j] > b || temp_lb[j] < b
                            temp_ub[j], temp_lb[j] = b, b
                            push!(I, j)
                        end
                    end
                end
            end
            #for implication bounds
            if imp_p_list[i] > 0
                record = imp_p[Int64(imp_p_list[i])]
                for rd in record
                    j = Int64(rd[1])
                    id = rd[2]
                    b = rd[3]
                    if id > 0
                        if temp_ub[j] > b
                            temp_ub[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
                        end
                    else
                        if temp_lb[j] < b
                            temp_lb[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
                        end
                    end
                end
            end
        else
            push!(J, i + binary_number)
            #change conflict bounds
            if temp_ub[bin_to_org[i]] > 0
                temp_ub[bin_to_org[i]] = 0
                push!(I, bin_to_org[i])
                (temp_lb[bin_to_org[i]] > 0) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
            end
            #for aggre
            if aggre_list[i] > 0
                record = aggre[Int64(aggre_list[i])] #here aggre includes 2 sparse vectors
                for rd in record
                    j = Int64(rd[1])
                    b = rd[3]
                    if temp_lb[j] > b || temp_ub[j] < b
                        return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]) #the last true implies we do not run this scenario
                    else
                        if temp_ub[j] > b || temp_lb[j] < b
                            temp_ub[j], temp_lb[j] = b, b
                            push!(I, j)
                        end
                    end
                end
            end
            #for implication bounds
            if imp_n_list[i] > 0
                record = imp_n[Int64(imp_n_list[i])]
                for rd in record
                    j = Int64(rd[1])
                    id = rd[2]
                    b = rd[3]
                    if id > 0
                        if temp_ub[j] > b
                            temp_ub[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
                        end
                    else
                        if temp_lb[j] < b
                            temp_lb[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
                        end
                    end
                end
            end
        end
    end
    
    #for two imp
    temp_imp_two = two_imp_list[J, J]
    K = findnz(temp_imp_two)[3]
    for k in K
        record = two_imp[Int64(k)]
        for rd in record
            i = Int64(rd[1])
            i_p = rd[2]
            b = rd[3]
            if i_p > 0
                if temp_ub[i] > b
                    temp_ub[i] = b
                    push!(I, i)
                    (temp_lb[i] > temp_ub[i]) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
                end
            else
                if temp_lb[i] < b
                    temp_lb[i] = b
                    push!(I, i)
                    (temp_lb[i] > temp_ub[i]) && (return temp_ub, temp_lb, spzeros(size(A)[1]), true, spzeros(size(A)[2]))
                end
            end
        end
    end
    
    #collect constraint numbers
    con_impact = spzeros(size(A)[1])
    var_impact = spzeros(size(A)[2])
    for i in I
        con_impact += A[:, i]
        var_impact[i] += 1
    end
    
    return temp_ub, temp_lb, con_impact, false, var_impact
end