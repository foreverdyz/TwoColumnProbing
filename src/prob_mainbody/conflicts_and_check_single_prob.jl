#conflicts_and_check_single_prob.jl

using SparseArrays

include("propagate_conflicts.jl")

"""
    conflicts_and_check_single_prob!()
Conduct conflict propagation and check whether fixing probing variable is feasible.
For instance, fixing x = 1 may make the problem to be infeasible. 
"""
function conflicts_and_check_single_prob!(
        x::Int64, y::Int64,
        var_ub::AbstractVector, var_lb::AbstractVector, bin_to_org::Dict{Int64, Int64},
        has_conflict_prop::AbstractVector, conflict_prop::AbstractVector, conflict_prop_count::Int64,
        clq_table::Vector{Vector{Int64}}, clq_p_n::Vector{Vector{Int64}}, S::Vector{Vector{Int64}}
    )    
    if var_ub[bin_to_org[y]] != var_lb[bin_to_org[y]] && var_ub[bin_to_org[x]] != var_lb[bin_to_org[x]]
        if has_conflict_prop[x] > 0
            x_cp, x_cn = conflict_prop[Int64(has_conflict_prop[x])]
            xp_id, xn_id = false, false #two falses implies x can be p and n
        else
            x_cp, xp_id = propagate_conflicts(x, true, clq_table, clq_p_n, S, length(bin_to_org))
            x_cn, xn_id = propagate_conflicts(x, false, clq_table, clq_p_n, S, length(bin_to_org))
            (xp_id && xn_id) && (println("problem infeasible"))
            conflict_prop_count += 1
            has_conflict_prop[x] = conflict_prop_count
            push!(conflict_prop, [x_cp, x_cn])
        end
        if has_conflict_prop[y] > 0
            y_cp, y_cn = conflict_prop[Int64(has_conflict_prop[y])]
            yp_id, yn_id = false, false
        else
            y_cp, yp_id = propagate_conflicts(y, true, clq_table, clq_p_n, S, length(bin_to_org))
            y_cn, yn_id = propagate_conflicts(y, false, clq_table, clq_p_n, S, length(bin_to_org))
            (yp_id && yn_id) && (println("problem infeasible"))
            conflict_prop_count += 1
            has_conflict_prop[y] = conflict_prop_count
            push!(conflict_prop, [y_cp, y_cn])
        end
        #get propagated conflicts
        return xp_id, xn_id, yp_id, yn_id, conflict_prop_count
    end
    
    if var_ub[bin_to_org[x]] == var_lb[bin_to_org[x]]
        if var_ub[bin_to_org[x]] > 0
            xp_id, xn_id = false, true # x can be 1 but cannot be 0
            if has_conflict_prop[x] > 0
                x_cp, _ = conflict_prop[Int64(has_conflict_prop[x])]
            else
                x_cp, xp_id = propagate_conflicts(x, true, clq_table, clq_p_n, S, length(bin_to_org))
                conflict_prop_count += 1
                has_conflict_prop[x] = conflict_prop_count
                push!(conflict_prop, [x_cp, spzeros(length(bin_to_org))])
            end
            (xp_id && xn_id) && (println("problem infeasible"))
        else
            xp_id, xn_id = true, false
            if has_conflict_prop[x] > 0
                _, x_cn = conflict_prop[Int64(has_conflict_prop[x])]
            else
                x_cn, xn_id = propagate_conflicts(x, false, clq_table, clq_p_n, S, length(bin_to_org))
                conflict_prop_count += 1
                has_conflict_prop[x] = conflict_prop_count
                push!(conflict_prop, [spzeros(length(bin_to_org)), x_cn])
            end
            (xp_id && xn_id) && (println("problem infeasible"))
        end
        if has_conflict_prop[y] > 0
            y_cp, y_cn = conflict_prop[Int64(has_conflict_prop[y])]
            yp_id, yn_id = false, false
        else
            y_cp, yp_id = propagate_conflicts(y, true, clq_table, clq_p_n, S, length(bin_to_org))
            y_cn, yn_id = propagate_conflicts(y, false, clq_table, clq_p_n, S, length(bin_to_org))
            (yp_id && yn_id) && (println("problem infeasible"))
            conflict_prop_count += 1
            has_conflict_prop[y] = conflict_prop_count
            push!(conflict_prop, [y_cp, y_cn])
        end
        #get propagated conflicts
        return xp_id, xn_id, yp_id, yn_id, conflict_prop_count
    elseif var_ub[bin_to_org[y]] == var_lb[bin_to_org[y]]
        if var_ub[bin_to_org[y]] > 0
            yp_id, yn_id = false, true
            if has_conflict_prop[y] > 0
                y_cp, _ = conflict_prop[Int64(has_conflict_prop[y])]
            else
                y_cp, yp_id = propagate_conflicts(y, true, clq_table, clq_p_n, S, length(bin_to_org))
                conflict_prop_count += 1
                has_conflict_prop[y] = conflict_prop_count
                push!(conflict_prop, [y_cp, spzeros(length(bin_to_org))])
            end
            (yp_id && yn_id) && (println("problem infeasible"))
        else
            yp_id, yn_id = true, false
            if has_conflict_prop[y] > 0
                _, y_cn = conflict_prop[Int64(has_conflict_prop[y])]
            else
                y_cn, yp_id = propagate_conflicts(y, false, clq_table, clq_p_n, S, length(bin_to_org))
                conflict_prop_count += 1
                has_conflict_prop[y] = conflict_prop_count
                push!(conflict_prop, [spzeros(length(bin_to_org)), y_cn])
            end
            (yp_id && yn_id) && (println("problem infeasible"))
        end
        if has_conflict_prop[x] > 0
            x_cp, x_cn = conflict_prop[Int64(has_conflict_prop[x])]
            xp_id, xn_id = false, false
        else
            x_cp, xp_id = propagate_conflicts(x, true, clq_table, clq_p_n, S, length(bin_to_org))
            x_cn, xn_id = propagate_conflicts(x, false, clq_table, clq_p_n, S, length(bin_to_org))
            (xp_id && xn_id) && (println("problem infeasible"))
            conflict_prop_count += 1
            has_conflict_prop[x] = conflict_prop_count
            push!(conflict_prop, [x_cp, x_cn])
        end
        #get propagated conflicts
        return xp_id, xn_id, yp_id, yn_id, conflict_prop_count
    end
end