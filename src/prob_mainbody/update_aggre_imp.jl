#update_aggre_imp.jl

"""
    update_aggre_imp!()
Update aggration form and implication graph based on probing detecting results.
"""
function update_aggre_imp!(
        z ,aggre_record, aggre, aggre_list, aggre_count,
        imp_p_ub, imp_p_lb, imp_p, imp_p_list, imp_p_count,
        imp_n_ub, imp_n_lb, imp_n, imp_n_list, imp_n_count,
    )
    if length(aggre_record) > 0
        if aggre_list[z] == 0
            aggre_count += 1
            aggre_list[z] = aggre_count
            push!(aggre, aggre_record)
        else
            j = Int64(aggre_list[z])
            append!!(aggre[j], aggre_record)
        end
    end
    if length(imp_p_ub) > 0
        if imp_p_list[z] == 0
            imp_p_count += 1
            imp_p_list[z] = imp_p_count
            push!(imp_p, imp_p_ub)
        else
            j = Int64(imp_p_list[z])
            append!!(imp_p[j], imp_p_ub)
        end
    end
    if length(imp_n_ub) > 0
        if imp_n_list[z] == 0
            imp_n_count += 1
            imp_n_list[z] = imp_n_count
            push!(imp_n, imp_n_ub)
        else
            j = Int64(imp_n_list[z])
            append!!(imp_n[j], imp_n_ub)
        end
    end
    if length(imp_p_lb) > 0
        if imp_p_list[z] == 0
            imp_p_count += 1
            imp_p_list[z] = imp_p_count
            push!(imp_p, imp_p_lb)
        else
            j = Int64(imp_p_list[z])
            append!!(imp_p[j], imp_p_lb)
        end
    end
    if length(imp_n_lb) > 0
        if imp_n_list[z] == 0
            imp_n_count += 1
            imp_n_list[z] = imp_n_count
            push!(imp_n, imp_n_lb)                        
        else
            j = Int64(imp_n_list[z])
            append!!(imp_n[j], imp_n_lb)
        end
    end
    return aggre_count, imp_p_count, imp_n_count
end