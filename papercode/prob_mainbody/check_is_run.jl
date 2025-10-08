#check_is_run.jl

using SparseArrays

"""
    check_is_run()
Check whether we start probing the two variables:
1, both two variables have been probed;
2, there are aggregation between two variables.
"""
function check_is_run(
        x::Int64, y::Int64, has_prob::AbstractVector, 
        var_ub::AbstractVector, var_lb::AbstractVector, bin_to_org::Dict{Int64, Int64},
        aggre_list::AbstractVector, is_aggred_list::AbstractVector, aggre::AbstractVector
    )
    #is_iter = min(has_prob[x] + is_aggred_list[x], 1) + min(has_prob[y] + is_aggred_list[y], 1)
    #if is_iter < 2 && var_ub[bin_to_org[y]] == var_lb[bin_to_org[y]] && var_ub[bin_to_org[x]] == var_lb[bin_to_org[x]]
    #    is_iter = 2
    #end
    
    is_iter = has_prob[x] + has_prob[y]
    if is_iter < 2 && var_ub[bin_to_org[y]] == var_lb[bin_to_org[y]] && var_ub[bin_to_org[x]] == var_lb[bin_to_org[x]]
        is_iter = 2
    end
    if is_iter < 2 && aggre_list[x] > 0
        record = aggre[Int64(aggre_list[x])]
        for rd in record
            if rd[1] == y
                is_iter = 2
                break
            end
        end
    end
    if is_iter < 2 && aggre_list[y] > 0
        record = aggre[Int64(aggre_list[y])]
        for rd in record
            if rd[1] == x
                is_iter = 2
                break
            end
        end
    end
    
    return is_iter
end