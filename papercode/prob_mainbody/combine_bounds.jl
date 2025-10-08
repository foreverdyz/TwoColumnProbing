#combine_bounds.jl

using SparseArrays


"""
    combine_bounds()
"""
function combine_bounds(
        x_ub::AbstractVector, x_lb::AbstractVector, 
        y_ub::AbstractVector, y_lb::AbstractVector, 
        x_cons::AbstractVector, y_cons::AbstractVector
    )
    impact_cons = findnz(x_cons + y_cons)[1]
    ub_temp = min.(x_ub, y_ub)
    lb_temp = max.(x_lb, y_lb)
    for i in 1:length(y_ub)
        (ub_temp[i] < lb_temp[i]) && (return ub_temp, lb_temp, impact_cons, false)
    end
    return ub_temp, lb_temp, impact_cons, true
end