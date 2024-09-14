#get_var_bounds.jl

using SparseArrays

include("get_temp_bounds.jl")

"""
    get_var_bounds()
Conduct conflict graph (clique table here) and implication graph propagation to get
implied variable bounds based on fixing probing variables.
"""
function get_var_bounds(
        x::Int64, y::Int64, xp_id::Bool, xn_id::Bool, yp_id::Bool, yn_id::Bool, 
        A::AbstractSparseArray, var_ub::AbstractVector, var_lb::AbstractVector, 
        conflict_prop::AbstractVector, has_conflict_prop::AbstractVector, 
        aggre_list::AbstractVector, aggre::AbstractVector,
        imp_p_list::AbstractVector, imp_p::AbstractVector, imp_n_list::AbstractVector, imp_n::AbstractVector,
        two_imp_list::AbstractSparseArray, two_imp::AbstractVector, binary_number::Int64, bin_to_org::Dict{Int64, Int64}
    )
    if xp_id
        xp_ub, xp_lb = spzeros(length(var_ub)), spzeros(length(var_lb))
        xp_impact = spzeros(size(A)[1])
    else
        xp_ub, xp_lb, xp_impact, xp_id, _ = get_temp_bounds(
                conflict_prop[Int64(has_conflict_prop[x])][1], binary_number,
                A, var_ub, var_lb,
                aggre_list, aggre,
                imp_p_list, imp_p, imp_n_list, imp_n,
                two_imp_list, two_imp, bin_to_org
            )
    end
    
    if xn_id
        xn_ub, xn_lb = spzeros(length(var_ub)), spzeros(length(var_lb))
        xn_impact = spzeros(size(A)[1])
        #push!(impact_set, xn_impact)
    else
        xn_ub, xn_lb, xn_impact, xn_id, _ = get_temp_bounds(
                conflict_prop[Int64(has_conflict_prop[x])][2], binary_number,
                A, var_ub, var_lb,
                aggre_list, aggre,
                imp_p_list, imp_p, imp_n_list, imp_n,
                two_imp_list, two_imp, bin_to_org
            )
        #push!(impact_set, xn_impact)
    end
    
    if yp_id
        yp_ub, yp_lb = spzeros(length(var_ub)), spzeros(length(var_lb))
        yp_impact = spzeros(size(A)[1])
        #push!(impact_set, yp_impact)
    else
        yp_ub, yp_lb, yp_impact, yp_id, _ = get_temp_bounds(
                conflict_prop[Int64(has_conflict_prop[y])][1], binary_number,
                A, var_ub, var_lb,
                aggre_list, aggre,
                imp_p_list, imp_p, imp_n_list, imp_n,
                two_imp_list, two_imp, bin_to_org
            )
        #push!(impact_set, yp_impact)
    end
    
    if yn_id
        yn_ub, yn_lb = spzeros(length(var_ub)), spzeros(length(var_lb))
        yn_impact = spzeros(size(A)[1])
        #push!(impact_set, yn_impact)
    else
        yn_ub, yn_lb, yn_impact, yn_id, _ = get_temp_bounds(
                conflict_prop[Int64(has_conflict_prop[y])][2], binary_number,
                A, var_ub, var_lb,
                aggre_list, aggre,
                imp_p_list, imp_p, imp_n_list, imp_n,
                two_imp_list, two_imp, bin_to_org
            )
        #push!(impact_set, yn_impact)
    end
    
    return xp_id, xn_id, yp_id, yn_id, xp_ub, xp_lb, xn_ub, xn_lb, yp_ub, yp_lb, yn_ub, yn_lb, xp_impact, xn_impact, yp_impact, yn_impact
end