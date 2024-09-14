#one_round_prob.jl

include("bounds_stren.jl")
include("combine_bounds.jl")

using SparseArrays

"""
    one_round_prob()
One iteration of two-column probing.
"""
function one_round_prob(
        is_xp::Bool, is_yp::Bool,
        xp_id::Bool, xn_id::Bool, yp_id::Bool, yn_id::Bool,
        x_ub::AbstractVector, x_lb::AbstractVector, y_ub::AbstractVector, y_lb::AbstractVector,
        x_cons::AbstractVector, y_cons::AbstractVector,
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, 
        con_ub::AbstractVector, con_lb::AbstractVector, var_type::AbstractVector
    )
    #check whether possible
    if is_xp
        #the last false implies we do not run this scenario
        (xp_id) && (return y_ub, y_lb, spzeros(length(y_ub)), spzeros(length(y_ub)), false)
    else
        #the last false implies we do not run this scenario
        (xn_id) && (return y_ub, y_lb, spzeros(length(y_ub)), spzeros(length(y_ub)), false)
    end
    
    if is_yp
        #the last false implies we do not run this scenario
        (yp_id) && (return x_ub, x_lb, spzeros(length(y_ub)), spzeros(length(y_ub)), false) 
    else
        #the last false implies we do not run this scenario
        (yn_id) && (return x_ub, x_lb, spzeros(length(y_ub)), spzeros(length(y_ub)), false)
    end
    
    #now start running
    #first, combine bounds
    ub_temp, lb_temp, impact_cons, run_prob_id = combine_bounds(x_ub, x_lb, y_ub, y_lb, x_cons, y_cons)
    if run_prob_id
        #second, run domain propagation
        return bounds_stren(impact_cons, con_set, con_coef, con_ub, con_lb, ub_temp, lb_temp, var_type)
    else
        #the last false implies we do not run this scenario
        return ub_temp, lb_temp, spzeros(length(ub_temp)), spzeros(length(lb_temp)), false
    end    
end