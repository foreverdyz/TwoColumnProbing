#update_conflict.jl

"""
    update_conflict()
Update conflict record, which is built in conflict graph (clique table here) propagation.
"""
function update_conflict(
        has_conflict_prop::AbstractVector, conflict_prop::AbstractVector,
        conflict_prop_count::Int64, new_conflict::AbstractVector,
        clq_table::Vector{Vector{Int64}}, clq_p_n::Vector{Vector{Int64}},
        S::Vector{Vector{Int64}}, binary_number::Int64
    )
    for rd in new_conflict
        x, y, x_id, y_id = Int64(rd[1]), Int64(rd[2]), rd[3], rd[4];
        if has_conflict_prop[x] <= 0
            x_cp, xp_id = propagate_conflicts(x, true, clq_table, clq_p_n, S, binary_number)
            x_cn, xn_id = propagate_conflicts(x, false, clq_table, clq_p_n, S, binary_number)
            (xp_id && xn_id) && (println("problem infeasible"))
            conflict_prop_count += 1
            has_conflict_prop[x] = conflict_prop_count
            push!(conflict_prop, [x_cp, x_cn])
        end
        if has_conflict_prop[y] <= 0
            y_cp, yp_id = propagate_conflicts(y, true, clq_table, clq_p_n, S, binary_number)
            y_cn, yn_id = propagate_conflicts(y, false, clq_table, clq_p_n, S, binary_number)
            (yp_id && yn_id) && (println("problem infeasible"))
            conflict_prop_count += 1
            has_conflict_prop[y] = conflict_prop_count
            push!(conflict_prop, [y_cp, y_cn])
        end
        if x_id > 0
            if y_id > 0
                conflict_prop[Int64(has_conflict_prop[x])][1] += conflict_prop[Int64(has_conflict_prop[y])][2]
                conflict_prop[Int64(has_conflict_prop[y])][1] += conflict_prop[Int64(has_conflict_prop[x])][2]
            else
                conflict_prop[Int64(has_conflict_prop[x])][1] += conflict_prop[Int64(has_conflict_prop[y])][1]
                conflict_prop[Int64(has_conflict_prop[y])][2] += conflict_prop[Int64(has_conflict_prop[x])][2]
            end
        else
            if y_id > 0
                conflict_prop[Int64(has_conflict_prop[x])][2] += conflict_prop[Int64(has_conflict_prop[y])][2]
                conflict_prop[Int64(has_conflict_prop[y])][1] += conflict_prop[Int64(has_conflict_prop[x])][1]
            else
                conflict_prop[Int64(has_conflict_prop[x])][2] += conflict_prop[Int64(has_conflict_prop[y])][1]
                conflict_prop[Int64(has_conflict_prop[y])][2] += conflict_prop[Int64(has_conflict_prop[x])][1]
            end
        end
    end
    return has_conflict_prop, conflict_prop, conflict_prop_count
end