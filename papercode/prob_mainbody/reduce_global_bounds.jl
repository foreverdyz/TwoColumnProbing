#reduce_global_bounds.jl

using Base.Threads

"""
    reduce_global_bounds()
After parallel two-column probing, we combine all bounds from different threads.
"""
function reduce_global_bounds(ub_g::AbstractVector, lb_g::AbstractVector, threadnum::Int64)
    let
        width = ceil(Int64, length(ub_g[1])/threadnum)
        @threads for id in 1:threadnum
            st = (id - 1) * width + 1
            en = min(id*width, length(ub_g[id]))
            for i in st:en
                ub_g[1][i] = minimum([ub_g[j][i] for j in 1:threadnum])
                lb_g[1][i] = maximum([lb_g[j][i] for j in 1:threadnum])
            end
        end
        return ub_g[1], lb_g[1]
    end
end