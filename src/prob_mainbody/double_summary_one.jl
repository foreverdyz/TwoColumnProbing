#double_summary_one.jl

using SparseArrays

"""
    double_summary_one()
Detect two-column implications.
"""
function double_summary_one(
        x::Int64, y::Int64, t_ub::AbstractVector, t_lb::AbstractVector, t_update_ub::AbstractVector, t_update_lb::AbstractVector,
        pc_ub::AbstractVector, pc_lb::AbstractVector, nc_ub::AbstractVector, nc_lb::AbstractVector
    )
    imp_rd = []
    #for ub
    for i in findnz(t_update_ub)[1]
        if t_ub[i] < pc_ub[i] && t_ub[i] < nc_ub[i]
            push!(imp_rd, [i, 1, t_ub[i]])
        end
    end
    for i in findnz(t_update_lb)[1]
        if t_lb[i] > pc_lb[i] && t_lb[i] > nc_lb[i]
            push!(imp_rd, [i, 0, t_lb[i]])
        end
    end
    return imp_rd
end