#presolve.jl

include("prepresolves/model_info.jl");
include("prepresolves/simple_presolve.jl");
include("prepresolves/classify_constraints.jl");
include("buildmatrices/clique_table.jl");
include("buildmatrices/coupling_matrix.jl");
include("prob_prepare/prob_order.jl");
include("prob_prepare/assign_variables.jl");
include("prob_prepare/prob_order_sub.jl");
include("prob_mainbody/probing.jl");
include("prob_mainbody/reduce_global_bounds.jl");
include("prob_mainbody/imp_analyze_global.jl");
include("rebuild_model/build_model_serial.jl");
include("rebuild_model/build_model_parallel.jl");

"""
    presolve()
Run our method to get a reduced MIP.
Intput:
filename: the name of .mps file
size_limit: the maximum size of constraints we considered in building coupling matrix
work_limit: the maximum number of elements processed in building coupling matrix
cand_number: the number of candidates we considered in probing order
max_probe_number: the maximum number of pairs of variables that we may probe (the upper bound of number of probing iterations)
time_limit: the time limit for our presolving
threadnum: the number of threads we used
path: the dir saving the original model
Output:
indicator: -1: time limit issue; 0: no impact model; 1: presolved model
time_used: runtime for presolving
model_reduce: reduced model
pair_number: number of pairs of variables have been probed
"""
function presolve(
        filename::String, size_limit::Int64, work_limit::Int64, 
        cand_number::Int64, max_probe_number::Int64, 
        time_limit::Real, threadnum::Int64, path::String
    )
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, 
    var_ub, var_lb, 
    var_name, var_type, 
    obj_coef, obj_constant, is_min = model_info(path*"/"*filename*".mps");
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, 
    var_ub, var_lb, var_type, 
    org_to_bin, bin_to_org = simple_presolve(
        con_matrix, con_set, con_coef, 
        con_ub, con_lb, var_ub, var_lb, var_type
    );
    I, S, B, J = classify_constraints(con_set, con_coef, con_ub, con_lb, var_type, org_to_bin);
    (length(I) < 1) && (return 0, 0, 0, 0)
    t = @timed begin
        start_time = time();
        clq_table, clq_p_n, clq_length = clique_table(S, length(bin_to_org));
        (time() - start_time > time_limit) && (return -1, 0, 0, 0)
        c_mt = coupling_matrix(B, length(org_to_bin), size_limit, work_limit);
        (time() - start_time > time_limit) && (return -1, 0, 0, 0)
        
        if threadnum > 1 #parallel
            assign_list, has_assigned = assign_variables(S, length(bin_to_org), threadnum);
            (time() - start_time > time_limit) && (return -1, -1, 0, 0)
            #initialize lists for different threads
            ub_g, lb_g = [var_ub for _ in 1:threadnum], [var_lb for _ in 1:threadnum]
            new_conflict_g = [[] for _ in 1:threadnum]
            aggre_list_g, aggre_g = [spzeros(length(bin_to_org)) for _ in 1:threadnum], [[] for _ in 1:threadnum]
            imp_p_list_g, imp_p_g = [spzeros(length(bin_to_org)) for _ in 1:threadnum], [[] for _ in 1:threadnum]
            imp_n_list_g, imp_n_g = [spzeros(length(bin_to_org)) for _ in 1:threadnum], [[] for _ in 1:threadnum]
            two_imp_list_g = [spzeros(2*length(bin_to_org), 2*length(bin_to_org)) for _ in 1:threadnum]
            two_imp_g = [[] for _ in 1:threadnum]
            has_conflict_prop_g, conflict_prop_g = [spzeros(length(bin_to_org)) for _ in 1:threadnum], [[] for _ in 1:threadnum]
            has_prob_g , conflict_prop_count_g= [spzeros(length(bin_to_org)) for _ in 1:threadnum], zeros(Int64, threadnum)
            pair_number_g = zeros(Int64, threadnum)
            len = zeros(threadnum);
            #start parallel two-column probing
            @threads for id in 1:threadnum
                #determine probing order
                prob_1, prob_2 = prob_order_sub(
                    con_matrix, clq_table, clq_length, 
                    c_mt, length(org_to_bin), assign_list[id], 
                    threadnum, bin_to_org, cand_number, max_probe_number
                );
                (time() - start_time > time_limit) && (return -1, 0, 0, 0)
                ub_g[id], lb_g[id], new_conflict_g[id], 
                aggre_list_g[id], _, aggre_g[id], 
                imp_p_list_g[id], imp_p_g[id], imp_n_list_g[id], imp_n_g[id], 
                two_imp_list_g[id], two_imp_g[id], has_conflict_prop_g[id], conflict_prop_g[id], conflict_prop_count_g[id],
                has_prob_g[id], pair_number_g[id] = probing(
                    prob_1, prob_2,
                    ub_g[id], lb_g[id], copy(var_type), 
                    copy(con_matrix[I, :]), copy(con_set[I]), 
                    copy(con_coef[I]), copy(con_ub[I]), copy(con_lb[I]),
                    copy(clq_table), copy(clq_p_n), copy(S), copy(bin_to_org), copy(org_to_bin), 
                    start_time, time_limit
                );
            end
            var_ub, var_lb = reduce_global_bounds(ub_g, lb_g, threadnum);
            var_ub, var_lb, aggre_g, aggre_list_g = imp_analyze_global(
                var_ub, var_lb, clq_table, clq_p_n, S, 
                bin_to_org, org_to_bin, aggre_list_g, aggre_g, 
                imp_p_list_g, imp_p_g, imp_n_list_g, imp_n_g, 
                two_imp_list_g, two_imp_g, has_conflict_prop_g, conflict_prop_g, 
                conflict_prop_count_g, has_assigned, has_prob_g, threadnum, start_time
            );
            pair_number = sum(pair_number_g);
        else #serial
            prob_1, prob_2 = prob_order(
                con_matrix, clq_table, clq_length, c_mt, 
                length(org_to_bin), bin_to_org, cand_number, max_probe_number
            );
            (time() - start_time > time_limit) && (return -1, 0, 0, 0)
            new_ub, new_lb, new_conflict, 
            aggre_list, is_aggred_list, aggre, 
            imp_p_list, imp_p,  imp_n_list, imp_n, two_imp_list, two_imp,
            has_conflict_prop, conflict_prop, conflict_prop_count,
            has_prob, pair_number = probing(
                prob_1, prob_2,
                copy(var_ub), copy(var_lb), var_type, 
                con_matrix[I, :], con_set[I], 
                con_coef[I], con_ub[I], con_lb[I],
                clq_table, clq_p_n, S, bin_to_org, org_to_bin, 
                start_time, time_limit
            );
        end
    end
    (pair_number < 1) && (return 0, 0, 0, 0)
    time_used = t.time - t.gctime;
    #rebuild model (get reduced model)
    if threadnum > 1 #parallel
        model_reduce = build_model_parallel(
            filename, var_ub, var_lb,
            new_conflict_g, aggre_list_g, aggre_g,
            bin_to_org, length(bin_to_org),
            is_min, obj_coef, threadnum
        );
    else #serial
        model_reduce = build_model_serial(
            filename, var_ub, var_lb,
            new_conflict, aggre_list, aggre,
            bin_to_org, length(bin_to_org),
            is_min, obj_coef
        );
    end
    return 1, time_used, model_reduce, pair_number
end