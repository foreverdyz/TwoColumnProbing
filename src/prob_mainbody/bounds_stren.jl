#bounds_stren.jl

using SparseArrays

"""
    bounds_stren()
Conduct domain propataion based on fixed variables from conflict graph (clique table here) and implication graph propagation.
"""
function bounds_stren(
        impact_cons::Vector{Int64}, 
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, 
        con_ub::AbstractVector, con_lb::AbstractVector, 
        ub_temp::AbstractVector, lb_temp::AbstractVector, var_type::AbstractVector
    )
    update_ub, update_lb = spzeros(length(ub_temp)), spzeros(length(lb_temp))
    for j in impact_cons
        if con_ub[j] < Inf
            stren_id = _ub_stren!(var_type, con_set[j], con_coef[j], con_ub[j], ub_temp, lb_temp, update_ub, update_lb)
            (stren_id) && (return ub_temp, lb_temp, update_ub, update_lb, false)
        end
        if con_lb[j] > -Inf
            stren_id = _lb_stren!(var_type, con_set[j], con_coef[j], con_lb[j], ub_temp, lb_temp, update_ub, update_lb)
            (stren_id) && (return ub_temp, lb_temp, update_ub, update_lb, false)
        end
    end
    return ub_temp, lb_temp, update_ub, update_lb, true 
end

function _ub_stren!(
        var_type::AbstractVector,
        con::Vector{Int64}, coef::AbstractVector, b::Real,
        var_ub::AbstractVector, var_lb::AbstractVector,
        update_ub::AbstractVector, update_lb::AbstractVector
    )
    #delta is the lower bound of LHS
    delta_list = zeros(length(con))
    for (i, j) in enumerate(con)
        #calculate lower bound for each variables
        if coef[i] > 0
            delta_list[i] = coef[i]*var_lb[j]
        else
            delta_list[i] = coef[i]*var_ub[j]
        end
    end
    delta = sum(delta_list)
    (delta > b + 0.0001) && (return true) #the last true means there is a conflict
    if delta != -Inf #avoid (-)Inf in variable bounds
        for (i, j) in enumerate(con)
            #the new potential bound
            local x = (b - delta + delta_list[i])/coef[i]
            if coef[i] > 0
                if var_ub[j] > x + 0.00001
                    (var_type[j] > 0) ? (var_ub[j] = floor(Int64, x)) : (var_ub[j] = x)
                    update_ub[j] = 1
                    #the last true means there is a conflict
                    (var_ub[j] < var_lb[j]) && (return true) 
                end
            else
                if var_lb[j] < x - 0.00001
                    (var_type[j] > 0) ? (var_lb[j] = ceil(Int64, x)) : (var_lb[j] = x)
                    update_lb[j] = 1
                    #the last true means there is a conflict
                    (var_ub[j] < var_lb[j]) && (return true)
                end
            end
        end
    end
    return false
end

function _lb_stren!(
        var_type::AbstractVector,
        con::Vector{Int64}, coef::AbstractVector, b::Real,
        var_ub::AbstractVector, var_lb::AbstractVector,
        update_ub::AbstractVector, update_lb::AbstractVector
    )
    #delta is the lower bound of LHS
    delta_list = zeros(length(con))
    for (i, j) in enumerate(con)
        #calculate lower bound for each variables
        if coef[i] > 0
            delta_list[i] = coef[i]*var_ub[j]
        else
            delta_list[i] = coef[i]*var_lb[j]
        end
    end
    delta = sum(delta_list)
    (delta < b - 0.0001) && (return true) #the last true means there is a conflict
    if delta != Inf #avoid (-)Inf in variable bounds
        for (i, j) in enumerate(con)
            #the new potential bound
            local x = (b - delta + delta_list[i])/coef[i]
            if coef[i] > 0
                if var_lb[j] < x + 0.00001
                    (var_type[j] > 0) ? (var_lb[j] = ceil(Int64, x)) : (var_lb[j] = x)
                    update_lb[j] = 1
                    #the last true means there is a conflict
                    (var_ub[j] < var_lb[j]) && (return true) 
                end
            else
                if var_ub[j] > x - 0.00001
                    (var_type[j] > 0) ? (var_ub[j] = floor(Int64, x)) : (var_ub[j] = x)
                    update_ub[j] = 1
                    #the last true means there is a conflict
                    (var_ub[j] < var_lb[j]) && (return true)
                end
            end
        end
    end
    return false
end