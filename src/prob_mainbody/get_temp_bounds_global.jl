#get_temp_bounds_global.jl

using SparseArrays

"""
    aggre[i]: [[j, i_p_j, i_n_j], ...]
    imp_p[i]: [[j, 1, j_ub] or [j, 0, j_lb], ...]
"""
function get_temp_bounds_global(
        conflict::AbstractVector, binary_number::Int64,
        var_ub::AbstractVector, var_lb::AbstractVector,
        aggre_list_g::AbstractVector, aggre_g::AbstractVector,
        imp_p_list_g::AbstractVector, imp_p_g::AbstractVector, 
        imp_n_list_g::AbstractVector, imp_n_g::AbstractVector,
        two_imp_list_g::AbstractVector, two_imp_g::AbstractVector,
        bin_to_org::Dict{Int64, Int64}, has_assigned::AbstractVector, threadnum::Int64
    )
    temp_ub, temp_lb = copy(var_ub), copy(var_lb)
    #construct temporary bounds
    I = Int64[]
    J = [Int64[] for _ in 1:threadnum]
    for i in findnz(conflict)[1]
        local id = Int64(has_assigned[i])
        if conflict[i] > 0
            push!(J[id], i)
            #change conflict bounds
            if temp_lb[bin_to_org[i]] < 1
                temp_lb[bin_to_org[i]] = 1
                push!(I, bin_to_org[i])
                (temp_ub[bin_to_org[i]] < 1) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
            end
            #for aggre
            if aggre_list_g[id][i] > 0
                record = aggre_g[id][Int64(aggre_list_g[id][i])]
                for rd in record
                    j = Int64(rd[1])
                    b = rd[2]
                    if temp_lb[j] > b || temp_ub[j] < b
                        return temp_ub, temp_lb, true, spzeros(length(var_ub)) 
                        #the last true implies we do not run this scenario
                    else
                        if temp_ub[j] > b || temp_lb[j] < b
                            temp_ub[j], temp_lb[j] = b, b
                            push!(I, j)
                        end
                    end
                end
            end
            #for implication bounds
            if imp_p_list_g[id][i] > 0
                record = imp_p_g[id][Int64(imp_p_list_g[id][i])]
                for rd in record
                    j = Int64(rd[1])
                    id = rd[2]
                    b = rd[3]
                    if id > 0
                        if temp_ub[j] > b
                            temp_ub[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
                        end
                    else
                        if temp_lb[j] < b
                            temp_lb[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
                        end
                    end
                end
            end
        else
            push!(J[id], i + binary_number)
            #change conflict bounds
            if temp_ub[bin_to_org[i]] > 0
                temp_ub[bin_to_org[i]] = 0
                push!(I, bin_to_org[i])
                (temp_lb[bin_to_org[i]] > 0) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
            end
            #for aggre
            local aggre_pos = Int64(aggre_list_g[id][i])
            if aggre_pos > 0
                record = aggre_g[id][aggre_pos] #here aggre includes 2 sparse vectors
                for rd in record
                    j = Int64(rd[1])
                    b = rd[3]
                    if temp_lb[j] > b || temp_ub[j] < b
                        return temp_ub, temp_lb, true, spzeros(length(var_ub))
                        #the last true implies we do not run this scenario
                    else
                        if temp_ub[j] > b || temp_lb[j] < b
                            temp_ub[j], temp_lb[j] = b, b
                            push!(I, j)
                        end
                    end
                end
            end
            #for implication bounds
            if imp_n_list_g[id][i] > 0
                record = imp_n_g[id][Int64(imp_n_list_g[id][i])]
                for rd in record
                    j = Int64(rd[1])
                    id = rd[2]
                    b = rd[3]
                    if id > 0
                        if temp_ub[j] > b
                            temp_ub[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
                        end
                    else
                        if temp_lb[j] < b
                            temp_lb[j] = b
                            push!(I, j)
                            (temp_lb[j] > temp_ub[j]) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
                        end
                    end
                end
            end
        end
    end
    
    #for two imp
    for id in 1:threadnum
        if length(J[id]) > 0
            temp_imp_two = two_imp_list_g[id][J[id], J[id]]
            K = findnz(temp_imp_two)[3]
            for k in K
                record = two_imp_g[id][Int64(k)]
                for rd in record
                    i = Int64(rd[1])
                    i_p = rd[2]
                    b = rd[3]
                    if i_p > 0
                        if temp_ub[i] > b
                            temp_ub[i] = b
                            push!(I, i)
                            (temp_lb[i] > temp_ub[i]) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
                        end
                    else
                        if temp_lb[i] < b
                            temp_lb[i] = b
                            push!(I, i)
                            (temp_lb[i] > temp_ub[i]) && (return temp_ub, temp_lb, true, spzeros(length(var_ub)))
                        end
                    end
                end
            end
        end 
    end
   
    #collect constraint numbers
    var_impact = spzeros(length(var_ub))
    for i in I
        var_impact[i] += 1
    end
    return temp_ub, temp_lb, false, var_impact
end