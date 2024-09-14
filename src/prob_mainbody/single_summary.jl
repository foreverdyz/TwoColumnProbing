#single_summary.jl

"""
    single_summary()
Dectect aggregation and one-column implication based on probing results.
"""
function single_summary(
        x::Int64,
        x_p_ub::AbstractVector, x_p_lb::AbstractVector, x_n_ub::AbstractVector, x_n_lb::AbstractVector,
        x_p_ub_update::AbstractVector, x_p_lb_update::AbstractVector, x_n_ub_update::AbstractVector, x_n_lb_update::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector
    )
    aggre_update_list = x_p_ub_update + x_p_lb_update + x_n_ub_update + x_n_lb_update
    aggre_record = []
    imp_p_ub = []
    imp_p_lb = []
    imp_n_ub = []
    imp_n_lb = []
    for i in findnz(aggre_update_list)[1]
        if i != x
            if x_p_ub[i] == x_p_lb[i] && x_n_ub[i] == x_n_lb[i]
                push!(aggre_record, [i, x_p_ub[i], x_n_ub[i]])
            else
                if x_p_ub[i] < var_ub[i]
                    push!(imp_p_ub, [i, 1, x_p_ub[i]])
                end
                if x_n_ub[i] < var_ub[i]
                    push!(imp_n_ub, [i, 1, x_n_ub[i]])
                end
                if x_p_lb[i] > var_lb[i]
                    push!(imp_p_lb, [i, 0, x_p_lb[i]])
                end
                if x_n_lb[i] > var_lb[i]
                    push!(imp_n_lb, [i, 0, x_n_lb[i]])
                end
            end
        end
    end
    return aggre_record, imp_p_ub, imp_p_lb, imp_n_ub, imp_n_lb
end