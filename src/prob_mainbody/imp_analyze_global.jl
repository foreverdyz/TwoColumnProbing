#imp_analyze_global.jl

include("get_temp_bounds_global.jl")

"""
    imp_analyze_global()
Conduct implication graph analysis to reduce/conclude probing results from different threads.
"""
function imp_analyze_global(
        var_ub::AbstractVector, var_lb::AbstractVector,
        clq_table::Vector{Vector{Int64}}, clq_p_n::AbstractVector, S::Vector{Vector{Int64}}, 
        bin_to_org::Dict{Int64, Int64}, org_to_bin::Dict{Int64, Int64},
        aggre_list_g::AbstractVector, aggre_g::AbstractVector, 
        imp_p_list_g::AbstractVector, imp_p_g::AbstractVector, 
        imp_n_list_g::AbstractVector, imp_n_g::AbstractVector, 
        two_imp_list_g::AbstractVector, two_imp_g::AbstractVector, 
        has_conflict_prop_g::AbstractVector, conflict_prop_g::AbstractVector, 
        conflict_prop_count_g::Vector{Int64},
        has_assigned::AbstractVector, has_prob_g::AbstractVector, threadnum::Int64, start_time::Real
    )
    inner_start_time = time();
    binary_number = length(bin_to_org);
    local ub_g, lb_g = [var_ub for _ in 1:threadnum], [var_lb for _ in 1:threadnum];
    for id in 1:threadnum
        for x in findnz(has_prob_g[id])[1]
            (time() - inner_start_time > 15) && (break)
            (time() - start_time > 75) && (break)
            xp_ub, xp_lb, xp_id, xp_impact = get_temp_bounds_global(
                conflict_prop_g[id][Int64(has_conflict_prop_g[id][x])][1], length(bin_to_org),
                var_ub, var_lb, aggre_list_g, aggre_g,
                imp_p_list_g, imp_p_g, imp_n_list_g, imp_n_g,
                two_imp_list_g, two_imp_g, bin_to_org, has_assigned, threadnum
            )
            xn_ub, xn_lb, xn_id, xn_impact = get_temp_bounds_global(
                conflict_prop_g[id][Int64(has_conflict_prop_g[id][x])][2], length(bin_to_org),
                var_ub, var_lb, aggre_list_g, aggre_g,
                imp_p_list_g, imp_p_g, imp_n_list_g, imp_n_g,
                two_imp_list_g, two_imp_g, bin_to_org, has_assigned, threadnum
            )
        
            if xp_id
                ub_g[id], lb_g[id] = xn_ub, xn_lb
            elseif xn_id
                ub_g[id], lb_g[id] = xp_ub, xp_lb
            else
                aggre_record = []
                impact_list = xp_impact + xn_impact
                for j in findnz(impact_list)[1]
                    if x != j
                        if xp_ub[j] == xp_lb[j] && xn_ub[j] == xn_lb[j]
                            push!(aggre_record, [j, xp_ub[j], xn_ub[j]])
                        else
                            if max(xp_ub[j], xn_ub[j]) < ub_g[id][j]
                                ub_g[id][j] = max(xp_ub[j], xn_ub[j])
                            end
                            if min(xp_lb[j], xn_lb[j]) > lb_g[id][j]
                                lb_g[id][j] = min(xp_lb[j], xn_lb[j])
                            end 
                        end
                    end
                end
                if length(aggre_record) > 0
                    aggre_count = length(findnz(aggre_list_g[id])[1])
                    if aggre_list_g[id][x] == 0
                        aggre_count += 1
                        aggre_list_g[id][x] = aggre_count
                        push!(aggre_g[id], aggre_record)
                    else
                        k = Int64(aggre_list_g[id][x])
                        append!!(aggre_g[id][k], aggre_record)
                    end
                end
            end
        end
    end
    var_ub, var_lb = reduce_global_bounds(ub_g, lb_g, threadnum);
    return var_ub, var_lb, aggre_g, aggre_list_g
end