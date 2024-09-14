#probing.jl

include("f_number_function.jl")
include("check_is_run.jl")
include("conflicts_and_check_single_prob.jl")
include("get_var_bounds.jl")
include("one_round_prob.jl")
include("single_summary.jl")
include("update_aggre_imp.jl")
include("double_summary_one.jl")
include("update_conflict.jl")

"""
    probing()
Two-colunm probing for pairs of binary variables (from prob_1 and prob_2 separately).
start_time is from previous steps, time_limit ia a customer parameter
Note that, we do not set max_probe_number here, because it is the length of prob_1 or prob_2
"""
function probing(
        prob_1::Vector{Int64}, prob_2::Vector{Int64},
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector,
        A::AbstractSparseArray, con_set::AbstractVector, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector,
        clq_table::Vector{Vector{Int64}}, clq_p_n::AbstractVector, S::Vector{Vector{Int64}},
        bin_to_org::Dict{Int64, Int64}, org_to_bin::Dict{Int64, Int64},
        start_time::Float64, time_limit::Real
    )
    #set max iteration
    max_iteration = length(prob_1);
    binary_number, var_num = length(bin_to_org), length(var_ub);
    
    #Initialize some sets
    #aggre is aggregated variables
    #aggre[i]: [j, i_p_j, i_n_j], i_p_j: i -> 1 then j value
    aggre_list, is_aggred_list, aggre_count, aggre = spzeros(binary_number), spzeros(binary_number), Int64(0), [];
    #imp is the implication graph
    #imp_p[i]: [j, 1, j_ub] or [j, 0, j_ub]
    imp_p_list, imp_p_count, imp_p = spzeros(binary_number), Int64(0), [];
    imp_n_list, imp_n_count, imp_n = spzeros(binary_number), Int64(0), [];
    #record whether this variable has been probbed
    has_prob = spzeros(binary_number);
    #record conflicts propagated
    has_conflict_prop, conflict_prop_count, conflict_prop = spzeros(binary_number), Int64(0), [];
    #new conflicts between binary variables
    #e.g. [i,j,1,1] implies i = 1 and j = 1 conflict
    new_conflict = [];
    #two columns imp
    two_imp_list, two_imp_count, two_imp = spzeros(2*binary_number, 2*binary_number), Int64(0), [];
    
    #mark efficiency
    eff = 0
    
    pair_number = 0
    
    for iter in 1:max_iteration
        if eff > 1000
            break
        end
        #get two variables to prob in this iteration
        x, y = prob_1[iter], prob_2[iter]
        #check whether we want to prob the two variables
        is_iter = check_is_run(
            x, y, has_prob, 
            var_ub, var_lb, bin_to_org,
            aggre_list, is_aggred_list, aggre
        )
        #start probing
        if is_iter < 2
            pair_number += 1
            eff = eff * 0.9
            #now we will prob the two variables separately
            has_prob[x] = 1
            has_prob[y] = 1
            #conflict propapate
            #xp_id = false implies x = 1 is feasible; otherwise x cannot be 1
            xp_id, xn_id, yp_id, yn_id, conflict_prop_count = conflicts_and_check_single_prob!(
                x, y, var_ub, var_lb, bin_to_org,
                has_conflict_prop, conflict_prop, conflict_prop_count,
                clq_table, clq_p_n, S
            )
            
            #get temp bounds (aggregation and implicationgraph propagate)
            xp_id, xn_id, yp_id, yn_id, 
            xp_ub, xp_lb, xn_ub, xn_lb, 
            yp_ub, yp_lb, yn_ub, yn_lb, 
            xp_impact, xn_impact, yp_impact, yn_impact = get_var_bounds(
                x, y, xp_id, xn_id, yp_id, yn_id, 
                A, var_ub, var_lb, 
                conflict_prop, has_conflict_prop, 
                aggre_list, aggre,
                imp_p_list, imp_p, imp_n_list, imp_n,
                two_imp_list, two_imp, binary_number, bin_to_org
            )
            
            #domain propagation
            pp_ub, pp_lb, pp_update_ub, pp_update_lb, pp_run_id = one_round_prob(
                true, true, xp_id, xn_id, yp_id, yn_id,
                xp_ub, xp_lb, yp_ub, yp_lb, xp_impact, yp_impact,
                con_set, con_coef, con_ub, con_lb, var_type
            )#1,1
            pn_ub, pn_lb, pn_update_ub, pn_update_lb, pn_run_id = one_round_prob(
                true, false, xp_id, xn_id, yp_id, yn_id,
                xp_ub, xp_lb, yn_ub, yn_lb, xp_impact, yn_impact,
                con_set, con_coef, con_ub, con_lb, var_type
            )#1,0
            np_ub, np_lb, np_update_ub, np_update_lb, np_run_id = one_round_prob(
                false, true, xp_id, xn_id, yp_id, yn_id,
                xn_ub, xn_lb, yp_ub, yp_lb, xn_impact, yp_impact,
                con_set, con_coef, con_ub, con_lb, var_type
            )#0,1
            nn_ub, nn_lb, nn_update_ub, nn_update_lb, nn_run_id = one_round_prob(
                false, false, xp_id, xn_id, yp_id, yn_id,
                xn_ub, xn_lb, yn_ub, yn_lb, xn_impact, yn_impact,
                con_set, con_coef, con_ub, con_lb, var_type
            )#0,0
            
            #collect and summarize results
            f_number = 0
            (!pp_run_id) && (f_number +=1)
            (!pn_run_id) && (f_number +=1)
            (!np_run_id) && (f_number +=1)
            (!nn_run_id) && (f_number +=1)
            if f_number == 4
                println("Problem Infeasible "*string(iter))
            elseif f_number == 3
                #here we fixed two binary variables
                eff = 0
                if pp_run_id
                    var_ub, var_lb = pp_ub, pp_lb
                elseif pn_run_id
                    var_ub, var_lb = pn_ub, pn_lb
                elseif np_run_id
                    var_ub, var_lb = np_ub, np_lb
                else
                    var_ub, var_lb = nn_ub, nn_lb
                end
            elseif f_number == 2
                eff = eff * 0.5;
                aggre_count, imp_p_count, imp_n_count, var_ub, var_lb = f_number_2!(
                    x, y,
                    var_ub, var_lb, var_type,
                    pp_run_id, pn_run_id, np_run_id, nn_run_id,
                    pp_ub, pp_lb, pn_ub, pn_lb, np_ub, np_lb, nn_ub, nn_lb,
                    pp_update_ub, pp_update_lb, pn_update_ub, pn_update_lb,
                    np_update_ub, np_update_lb, nn_update_ub, nn_update_lb,
                    aggre, aggre_list, is_aggred_list, aggre_count,
                    imp_p, imp_p_list, imp_p_count,
                    imp_n, imp_n_list, imp_n_count,
                    org_to_bin
                );                
            elseif f_number == 1
                eff = eff * 0.8;
                aggre_count, imp_p_count, imp_n_count, two_imp_count, var_ub, var_lb = f_number_1!(
                    x, y,
                    var_ub, var_lb, var_type,
                    pp_run_id, pn_run_id, np_run_id, nn_run_id,
                    pp_ub, pp_lb, pn_ub, pn_lb, np_ub, np_lb, nn_ub, nn_lb,
                    pp_update_ub, pp_update_lb, pn_update_ub, pn_update_lb,
                    np_update_ub, np_update_lb, nn_update_ub, nn_update_lb,
                    aggre, aggre_list, is_aggred_list, aggre_count,
                    imp_p, imp_p_list, imp_p_count,
                    imp_n, imp_n_list, imp_n_count,
                    two_imp_list, two_imp, two_imp_count,
                    new_conflict, binary_number, org_to_bin
                );
                #has_conflict_prop, conflict_prop, conflict_prop_count = update_conflict(
                #    has_conflict_prop, conflict_prop, 
                #    conflict_prop_count, [new_conflict[end]],
                #    clq_table, clq_p_n, S, length(bin_to_org)
                #);
            else
                aggre_count, imp_p_count, imp_n_count, two_imp_count, eff, var_ub, var_lb = f_number_0!(
                    x, y,
                    var_ub, var_lb, var_type,
                    pp_run_id, pn_run_id, np_run_id, nn_run_id,
                    pp_ub, pp_lb, pn_ub, pn_lb, np_ub, np_lb, nn_ub, nn_lb,
                    pp_update_ub, pp_update_lb, pn_update_ub, pn_update_lb,
                    np_update_ub, np_update_lb, nn_update_ub, nn_update_lb,
                    aggre, aggre_list, is_aggred_list, aggre_count,
                    imp_p, imp_p_list, imp_p_count,
                    imp_n, imp_n_list, imp_n_count,
                    two_imp_list, two_imp, two_imp_count,
                    eff, binary_number, org_to_bin
                );
            end
            eff += 110
        end
        if time() - start_time > time_limit
            break
        end
    end
    return var_ub, var_lb, new_conflict, aggre_list, is_aggred_list, aggre, imp_p_list, imp_p, imp_n_list, imp_n, two_imp_list, two_imp, has_conflict_prop, conflict_prop, conflict_prop_count, has_prob, pair_number
end