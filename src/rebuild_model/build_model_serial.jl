#build_model_serial.jl

using JuMP
using SparseArrays

function build_model_serial(
        filename, var_ub, var_lb,
        new_conflict, aggre_list, aggre,
        bin_to_org, binary_number,
        is_min, obj_coef
    )
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, _, _, 
    var_name, var_type, 
    obj_coef, obj_constant, is_min = model_info("mipdata/"*filename*".mps");
    con_num = length(con_ub);
    con_type = zeros(con_num);
    var_num = length(var_ub);
    
    model = Model()
    
    #build variables
    @variable(model, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])

    #set binary and integer variables
    for i in 1:var_num
        (var_type[i] == 2) && (set_binary(x[i]))
        (var_type[i] == 1) && (set_integer(x[i]))
    end
    
    for j in 1:con_num
        if con_ub[j] == con_lb[j]
            @constraint(model, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) <= con_ub[j])
            @constraint(model, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) >= con_lb[j])
        else
            if con_ub[j] < Inf
                @constraint(model, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) <= con_ub[j])
            end
            if con_lb[j] > -Inf
                @constraint(model, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) >= con_lb[j])
            end
        end    
    end
    
    if length(new_conflict) > 0
        for rd in new_conflict
            i = bin_to_org[Int64(rd[1])]
            j = bin_to_org[Int64(rd[2])]
            c_i = 2*(rd[3] - 0.5)
            c_j = 2*(rd[4] - 0.5)
            b = 1 - 2 + rd[3] + rd[4]
            @constraint(model, c_i*x[i] + c_j*x[j] <= b)
        end
    end
    
    if length(aggre) > 0
        for i in findnz(aggre_list)[1]
            record = aggre[Int64(aggre_list[i])]
            i_org = bin_to_org[i]
            for rd in record
                j = Int64(rd[1])
                p = rd[2]
                n = rd[3]
                @constraint(model, p*x[i_org] + n*(1-x[i_org]) - x[j] == 0)
            end
        end
    end
    #set objective function
    (is_min) ? (@objective(model, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))
    return model
end