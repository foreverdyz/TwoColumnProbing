#f_number_function.jl


"""
    Define 3 functions for there are 0, 1, 2 conflicts between two binary variables that we are probing.
"""

function f_number_0!(
        x::Int64, y::Int64,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector,
        pp_run_id::Bool, pn_run_id::Bool, np_run_id::Bool, nn_run_id::Bool,
        pp_ub::AbstractVector, pp_lb::AbstractVector, pn_ub::AbstractVector, pn_lb::AbstractVector,
        np_ub::AbstractVector, np_lb::AbstractVector, nn_ub::AbstractVector, nn_lb::AbstractVector,
        pp_update_ub::AbstractVector, pp_update_lb::AbstractVector, 
        pn_update_ub::AbstractVector, pn_update_lb::AbstractVector,
        np_update_ub::AbstractVector, np_update_lb::AbstractVector, 
        nn_update_ub::AbstractVector, nn_update_lb::AbstractVector,
        aggre::AbstractVector, aggre_list::AbstractVector, is_aggred_list::AbstractVector, aggre_count::Int64,
        imp_p::AbstractVector, imp_p_list::AbstractVector, imp_p_count::Int64,
        imp_n::AbstractVector, imp_n_list::AbstractVector, imp_n_count::Int64,
        two_imp_list::SparseMatrixCSC, two_imp::AbstractVector, two_imp_count::Int64,
        eff::Real, binary_number::Int64, org_to_bin::Dict{Int64, Int64}
    )
    var_ub, var_lb = min.(var_ub, max.(pp_ub, np_ub, pn_ub, nn_ub)), max.(var_lb, min.(pp_lb, np_lb, pn_lb, nn_lb))
    x_p_ub, x_p_lb, x_n_ub, x_n_lb = max.(pp_ub, pn_ub), min.(pp_lb, pn_lb), max.(np_ub, nn_ub), min.(np_lb, nn_lb)
    x_p_ub_list, x_p_lb_list = pn_update_ub + pp_update_ub, pn_update_lb + pp_update_lb
    x_n_ub_list, x_n_lb_list = np_update_ub + nn_update_ub, np_update_lb + nn_update_lb
    y_p_ub, y_p_lb, y_n_ub, y_n_lb = max.(pp_ub, np_ub), min.(pp_lb, np_lb), max.(pn_ub, nn_ub), min.(pn_lb, nn_lb)
    y_p_ub_list, y_p_lb_list = np_update_ub + pp_update_ub, np_update_lb + pp_update_lb
    y_n_ub_list, y_n_lb_list = pn_update_ub + nn_update_ub, pn_update_lb + nn_update_lb
         
    #update aggregation and imp graph
    aggre_record, imp_p_ub, imp_p_lb, imp_n_ub, imp_n_lb = single_summary(
            x, 
            x_p_ub, x_p_lb, x_n_ub, x_n_lb, 
            x_p_ub_list, x_p_lb_list, x_n_ub_list, x_n_lb_list,
            var_ub, var_lb
    )
    aggre_count, imp_p_count, imp_n_count = update_aggre_imp!(
        x, aggre_record, aggre, aggre_list, aggre_count,
        imp_p_ub, imp_p_lb, imp_p, imp_p_list, imp_p_count,
        imp_n_ub, imp_n_lb, imp_n, imp_n_list, imp_n_count,
    )
    for rd in aggre_record
        z = Int64(rd[1])
        if var_type[z] >= 2
            is_aggred_list[org_to_bin[z]] = 1
        end
    end
    if length(aggre_record) > 1
        eff = eff * 0.9
    #elseif length(imp_p_ub) + length(imp_p_lb) + length(imp_n_ub) + length(imp_n_lb)
    #    eff = eff * 0.95
    end
    
    aggre_record, imp_p_ub, imp_p_lb, imp_n_ub, imp_n_lb = single_summary(
        y, 
        y_p_ub, y_p_lb, y_n_ub, y_n_lb, 
        y_p_ub_list, x_p_lb_list, y_n_ub_list, y_n_lb_list,
        var_ub, var_lb
    )
    aggre_count, imp_p_count, imp_n_count = update_aggre_imp!(
        y, aggre_record, aggre, aggre_list, aggre_count,
        imp_p_ub, imp_p_lb, imp_p, imp_p_list, imp_p_count,
        imp_n_ub, imp_n_lb, imp_n, imp_n_list, imp_n_count,
    )
    for rd in aggre_record
        z = Int64(rd[1])
        if var_type[z] >= 2
            is_aggred_list[org_to_bin[z]] = 1
        end
    end
    if length(aggre_record) > 1
        eff = eff * 0.9
    #elseif length(imp_p_ub) + length(imp_p_lb) + length(imp_n_ub) + length(imp_n_lb)
    #    eff = eff * 0.95
    end
    
    #update two imp graph
    #for pp
    imp_rd = double_summary_one(x, y, pp_ub, pp_lb, pp_update_ub, pp_update_lb, pn_ub, pn_lb, np_ub, np_lb)
    if length(imp_rd) > 0
        if two_imp_list[x, y] <= 0
            two_imp_count += 1
            two_imp_list[x, y] = two_imp_count
            push!(two_imp, imp_rd)
        else
            append!!(two_imp[Int64(two_imp_list[x, y])], imp_rd)
        end
    end
    if length(imp_rd) > 1
        eff = eff * 0.95
    end
    #for pn
    imp_rd = double_summary_one(x, y, pn_ub, pn_lb, pn_update_ub, pn_update_lb, pp_ub, pp_lb, nn_ub, nn_lb)
    if length(imp_rd) > 0
        if two_imp_list[x, y + binary_number] <= 0
            two_imp_count += 1
            two_imp_list[x, y + binary_number] = two_imp_count
            push!(two_imp, imp_rd)
        else
            append!!(two_imp[Int64(two_imp_list[x, y + binary_number])], imp_rd)
        end
    end
    if length(imp_rd) > 1
        eff = eff * 0.95
    end
    #for np
    imp_rd = double_summary_one(x, y, np_ub, np_lb, np_update_ub, np_update_lb, pp_ub, pp_lb, nn_ub, nn_lb)
    if length(imp_rd) > 0
        if two_imp_list[x + binary_number, y] <= 0
            two_imp_count += 1
            two_imp_list[x + binary_number, y] = two_imp_count
            push!(two_imp, imp_rd)
        else
            append!!(two_imp[Int64(two_imp_list[x + binary_number, y])], imp_rd)
        end
    end
    if length(imp_rd) > 1
        eff = eff * 0.95
    end
    #for nn
    imp_rd = double_summary_one(x, y, nn_ub, nn_lb, nn_update_ub, nn_update_lb, pn_ub, pn_lb, np_ub, np_lb)
    if length(imp_rd) > 0
        if two_imp_list[x + binary_number, y + binary_number] <= 0
            two_imp_count += 1
            two_imp_list[x + binary_number, y + binary_number] = two_imp_count
            push!(two_imp, imp_rd)
        else
            append!!(two_imp[Int64(two_imp_list[x + binary_number, y + binary_number])], imp_rd)
        end
    end
    if length(imp_rd) > 1
        eff = eff * 0.95
    end
    
    return aggre_count, imp_p_count, imp_n_count, two_imp_count, eff, var_ub, var_lb
end

function f_number_1!(
        x::Int64, y::Int64,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector,
        pp_run_id::Bool, pn_run_id::Bool, np_run_id::Bool, nn_run_id::Bool,
        pp_ub::AbstractVector, pp_lb::AbstractVector, pn_ub::AbstractVector, pn_lb::AbstractVector,
        np_ub::AbstractVector, np_lb::AbstractVector, nn_ub::AbstractVector, nn_lb::AbstractVector,
        pp_update_ub::AbstractVector, pp_update_lb::AbstractVector, 
        pn_update_ub::AbstractVector, pn_update_lb::AbstractVector,
        np_update_ub::AbstractVector, np_update_lb::AbstractVector, 
        nn_update_ub::AbstractVector, nn_update_lb::AbstractVector,
        aggre::AbstractVector, aggre_list::AbstractVector, is_aggred_list::AbstractVector, aggre_count::Int64,
        imp_p::AbstractVector, imp_p_list::AbstractVector, imp_p_count::Int64,
        imp_n::AbstractVector, imp_n_list::AbstractVector, imp_n_count::Int64,
        two_imp_list::SparseMatrixCSC, two_imp::AbstractVector, two_imp_count::Int64,
        new_conflict::AbstractVector, binary_number::Int64,  org_to_bin::Dict{Int64, Int64}
    )
    #collect implied bounds for x and y and update two imp graph
    if !pp_run_id
        var_ub, var_lb = min.(var_ub, max.(np_ub, pn_ub, nn_ub)), max.(var_lb, min.(np_lb, pn_lb, nn_lb));
        x_p_ub, x_p_lb, x_n_ub, x_n_lb = pn_ub, pn_lb, max.(np_ub, nn_ub), min.(np_lb, nn_lb);
        x_p_ub_list, x_p_lb_list = pn_update_ub, pn_update_lb;
        x_n_ub_list, x_n_lb_list = np_update_ub + nn_update_ub, np_update_lb + nn_update_lb;
        y_p_ub, y_p_lb, y_n_ub, y_n_lb = np_ub, np_lb, max.(pn_ub, nn_ub), min.(pn_lb, nn_lb);
        y_p_ub_list, y_p_lb_list = np_update_ub, np_update_lb;
        y_n_ub_list, y_n_lb_list = pn_update_ub + nn_update_ub, pn_update_lb + nn_update_lb;
        #for nn
        imp_rd = double_summary_one(x, y, nn_ub, nn_lb, nn_update_ub, nn_update_lb, pn_ub, pn_lb, np_ub, np_lb)
        if length(imp_rd) > 0
            if two_imp_list[x + binary_number, y + binary_number] <= 0
                two_imp_count += 1
                two_imp_list[x + binary_number, y + binary_number] = two_imp_count
                push!(two_imp, imp_rd)
            else
                append!!(two_imp[Int64(two_imp_list[x + binary_number, y + binary_number])], imp_rd)
            end
        end
        push!(new_conflict, [x, y, 1, 1])
    elseif !pn_run_id
        var_ub, var_lb = min.(var_ub, max.(np_ub, pp_ub, nn_ub)), max.(var_lb, min.(np_lb, pp_lb, nn_lb))
        x_p_ub, x_p_lb, x_n_ub, x_n_lb = pp_ub, pp_lb, max.(np_ub, nn_ub), min.(np_lb, nn_lb)
        x_p_ub_list, x_p_lb_list = pp_update_ub, pp_update_lb
        x_n_ub_list, x_n_lb_list = np_update_ub + nn_update_ub, np_update_lb + nn_update_lb
        y_p_ub, y_p_lb, y_n_ub, y_n_lb = max.(pp_ub, np_ub), min.(pp_lb, np_lb), nn_ub, nn_lb
        y_p_ub_list, y_p_lb_list = np_update_ub + pp_update_ub, np_update_lb + pp_update_lb
        y_n_ub_list, y_n_lb_list = nn_update_ub, nn_update_lb
        #for np
        imp_rd = double_summary_one(x, y, np_ub, np_lb, np_update_ub, np_update_lb, pp_ub, pp_lb, nn_ub, nn_lb)
        if length(imp_rd) > 0
            if two_imp_list[x + binary_number, y] <= 0
                two_imp_count += 1
                two_imp_list[x + binary_number, y] = two_imp_count
                push!(two_imp, imp_rd)
            else
                append!!(two_imp[Int64(two_imp_list[x + binary_number, y])], imp_rd)
            end
        end
        push!(new_conflict, [x, y, 1, 0])
    elseif !np_run_id
        var_ub, var_lb = min.(var_ub, max.(pn_ub, pp_ub, nn_ub)), max.(var_lb, min.(pn_lb, pp_lb, nn_lb))
        x_p_ub, x_p_lb, x_n_ub, x_n_lb = max.(pp_ub, pn_ub), min.(pp_lb, pn_lb), nn_ub, nn_lb
        x_p_ub_list, x_p_lb_list = pp_update_ub + pn_update_ub, pp_update_lb + pn_update_lb
        x_n_ub_list, x_n_lb_list = nn_update_ub, nn_update_lb
        y_p_ub, y_p_lb, y_n_ub, y_n_lb = pp_ub, pp_lb, max.(pn_ub, nn_ub), min.(pn_lb, nn_lb)
        y_p_ub_list, y_p_lb_list = pp_update_ub, pp_update_lb
        y_n_ub_list, y_n_lb_list = pn_update_ub + nn_update_ub, pn_update_lb + nn_update_lb
        #for pn
        imp_rd = double_summary_one(x, y, pn_ub, pn_lb, pn_update_ub, pn_update_lb, pp_ub, pp_lb, nn_ub, nn_lb)
        if length(imp_rd) > 0
            if two_imp_list[x, y + binary_number] <= 0
                two_imp_count += 1
                two_imp_list[x, y + binary_number] = two_imp_count
                push!(two_imp, imp_rd)
            else
                append!!(two_imp[Int64(two_imp_list[x, y + binary_number])], imp_rd)
            end
        end
        push!(new_conflict, [x, y, 0, 1])
    else
        var_ub, var_lb= min.(var_ub, max.(pn_ub, pp_ub, np_ub)), max.(var_lb, min.(pn_lb, pp_lb, np_lb))
        x_p_ub, x_p_lb, x_n_ub, x_n_lb = max.(pp_ub, pn_ub), min.(pp_lb, pn_lb), np_ub, np_lb
        x_p_ub_list, x_p_lb_list = pp_update_ub + pn_update_ub, pp_update_lb + pn_update_lb
        x_n_ub_list, x_n_lb_list = np_update_ub, np_update_lb
        y_p_ub, y_p_lb, y_n_ub, y_n_lb = max.(pp_ub, np_ub), min.(pp_lb, np_lb), pn_ub, pn_lb
        y_p_ub_list, y_p_lb_list = pp_update_ub + np_update_ub, pp_update_lb + np_update_lb
        y_n_ub_list, y_n_lb_list = pn_update_ub, pn_update_lb
        #for pp
        imp_rd = double_summary_one(x, y, pp_ub, pp_lb, pp_update_ub, pp_update_lb, pn_ub, pn_lb, np_ub, np_lb)
        if length(imp_rd) > 0
            if two_imp_list[x, y] <= 0
                two_imp_count += 1
                two_imp_list[x, y] = two_imp_count
                push!(two_imp, imp_rd)
            else
                append!!(two_imp[Int64(two_imp_list[x, y])], imp_rd)
            end
        end
        push!(new_conflict, [x, y, 0, 0])
    end
    
    #update aggregate and imp graph for x
    aggre_record, imp_p_ub, imp_p_lb, imp_n_ub, imp_n_lb = single_summary(
            x, 
            x_p_ub, x_p_lb, x_n_ub, x_n_lb, 
            x_p_ub_list, x_p_lb_list, x_n_ub_list, x_n_lb_list,
            var_ub, var_lb
    )
    aggre_count, imp_p_count, imp_n_count = update_aggre_imp!(
        x, aggre_record, aggre, aggre_list, aggre_count,
        imp_p_ub, imp_p_lb, imp_p, imp_p_list, imp_p_count,
        imp_n_ub, imp_n_lb, imp_n, imp_n_list, imp_n_count,
    )
    for rd in aggre_record
        z = Int64(rd[1])
        if var_type[z] >= 2
            is_aggred_list[org_to_bin[z]] = 1
        end
    end
    #update aggregate and imp graph for y
    aggre_record, imp_p_ub, imp_p_lb, imp_n_ub, imp_n_lb = single_summary(
            y, 
            y_p_ub, y_p_lb, y_n_ub, y_n_lb, 
            y_p_ub_list, x_p_lb_list, y_n_ub_list, y_n_lb_list,
            var_ub, var_lb
    )
    aggre_count, imp_p_count, imp_n_count = update_aggre_imp!(
        y, aggre_record, aggre, aggre_list, aggre_count,
        imp_p_ub, imp_p_lb, imp_p, imp_p_list, imp_p_count,
        imp_n_ub, imp_n_lb, imp_n, imp_n_list, imp_n_count,
    )
    for rd in aggre_record
        z = Int64(rd[1])
        if var_type[z] >= 2
            is_aggred_list[org_to_bin[z]] = 1
        end
    end
    return aggre_count, imp_p_count, imp_n_count, two_imp_count, var_ub, var_lb
end

function f_number_2!(
        x::Int64, y::Int64,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector,
        pp_run_id::Bool, pn_run_id::Bool, np_run_id::Bool, nn_run_id::Bool,
        pp_ub::AbstractVector, pp_lb::AbstractVector, pn_ub::AbstractVector, pn_lb::AbstractVector,
        np_ub::AbstractVector, np_lb::AbstractVector, nn_ub::AbstractVector, nn_lb::AbstractVector,
        pp_update_ub::AbstractVector, pp_update_lb::AbstractVector, 
        pn_update_ub::AbstractVector, pn_update_lb::AbstractVector,
        np_update_ub::AbstractVector, np_update_lb::AbstractVector, 
        nn_update_ub::AbstractVector, nn_update_lb::AbstractVector,
        aggre::AbstractVector, aggre_list::AbstractVector, is_aggred_list::AbstractVector, aggre_count::Int64,
        imp_p::AbstractVector, imp_p_list::AbstractVector, imp_p_count::Int64,
        imp_n::AbstractVector, imp_n_list::AbstractVector, imp_n_count::Int64,
        org_to_bin::Dict{Int64, Int64}
    )
    z = 0
    var_num = length(var_ub)
    p_ub, p_lb, n_ub, n_lb = spzeros(var_num), spzeros(var_num), spzeros(var_num), spzeros(var_num);
    p_update_ub, p_update_lb = spzeros(var_num), spzeros(var_num);
    n_update_ub, n_update_lb = spzeros(var_num), spzeros(var_num);
    if pp_run_id && pn_run_id
        var_ub, var_lb = min.(var_ub, max.(pp_ub, pn_ub)), max.(var_lb, min.(pp_lb, pn_lb));
        z = y;
        p_ub, p_lb, n_ub, n_lb = pp_ub, pp_lb, pn_ub, pn_lb;
        p_update_ub, p_update_lb, n_update_ub, n_update_lb = pp_update_ub, pp_update_lb, pn_update_ub, pn_update_lb;
    elseif pp_run_id && np_run_id
        var_ub, var_lb = min.(var_ub, max.(pp_ub, np_ub)), max.(var_lb, min.(pp_lb, np_lb));
        z = x;
        p_ub, p_lb, n_ub, n_lb = pp_ub, pp_lb, np_ub, np_lb;
        p_update_ub, p_update_lb, n_update_ub, n_update_lb = pp_update_ub, pp_update_lb, np_update_ub, np_update_lb;
    elseif pp_run_id && nn_run_id
        var_ub, var_lb = min.(var_ub, max.(pp_ub, nn_ub)), max.(var_lb, min.(pp_lb, nn_lb));
        z = x;
        p_ub, p_lb, n_ub, n_lb = pp_ub, pp_lb, nn_ub, nn_lb;
        p_update_ub, p_update_lb, n_update_ub, n_update_lb = pp_update_ub, pp_update_lb, nn_update_ub, nn_update_lb;
    elseif pn_run_id && np_run_id
        var_ub, var_lb = min.(var_ub, max.(pn_ub, np_ub)), max.(var_lb, min.(pn_lb, np_lb));
        z = x;
        p_ub, p_lb, n_ub, n_lb = pn_ub, pn_lb, np_ub, np_lb;
        p_update_ub, p_update_lb, n_update_ub, n_update_lb = pn_update_ub, pn_update_lb, np_update_ub, np_update_lb;
    elseif pn_run_id && nn_run_id
        var_ub, var_lb = min.(var_ub, max.(pn_ub, nn_ub)), max.(var_lb, min.(pn_lb, nn_lb));
        z = x;
        p_ub, p_lb, n_ub, n_lb = pn_ub, pn_lb, nn_ub, nn_lb;
        p_update_ub, p_update_lb, n_update_ub, n_update_lb = pn_update_ub, pn_update_lb, nn_update_ub, nn_update_lb;
    else
        var_ub, var_lb = min.(var_ub, max.(np_ub, nn_ub)), max.(var_lb, min.(np_lb, nn_lb));
        z = y;
        p_ub, p_lb, n_ub, n_lb = np_ub, np_lb, nn_ub, nn_lb;
        p_update_ub, p_update_lb, n_update_ub, n_update_lb = np_update_ub, np_update_lb, nn_update_ub, nn_update_lb;
    end
    aggre_record, imp_p_ub, imp_p_lb, imp_n_ub, imp_n_lb = single_summary(
            z, 
            p_ub, p_lb, n_ub, n_lb, 
            p_update_ub, p_update_lb, n_update_ub, n_update_lb,
            var_ub, var_lb
    );
    aggre_count, imp_p_count, imp_n_count = update_aggre_imp!(
            z, aggre_record, aggre, aggre_list, aggre_count,
            imp_p_ub, imp_p_lb, imp_p, imp_p_list, imp_p_count,
            imp_n_ub, imp_n_lb, imp_n, imp_n_list, imp_n_count,
    );
    for rd in aggre_record
        z = Int64(rd[1])
        if var_type[z] >= 2
            is_aggred_list[org_to_bin[z]] = 1
        end
    end
    return aggre_count, imp_p_count, imp_n_count, var_ub, var_lb
end