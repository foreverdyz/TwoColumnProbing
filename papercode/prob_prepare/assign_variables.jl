#assign_variables.jl

using SparseArrays

#here we split varaibles: one clique in one thread
function assign_variables(S::Vector{Vector{Int64}}, binary_number::Int64, threadnum::Int64)
    #build a sparse list to cache whether a binary varaible has been assigned
    has_assigned = spzeros(binary_number)
    #pos is the next to-be-assigned thread
    pos = 1
    #assigned_list cache index in each thread
    assigned_list = [Int64[] for _ in 1:threadnum]
    assigned_len = Int64[length(assigned_list[i]) for i in 1:threadnum]
    #total_assigned denotes how many varaibles have been assigned
    total_assigned = 0
    #start assignment based on cliques
    for (i, s) in enumerate(S)
        for x in s
            (x > binary_number) && (x -= binary_number)
            if has_assigned[x] < 1
                assigned_len[pos] += 1
                has_assigned[x] = pos
                push!(assigned_list[pos], x)
                total_assigned += 1
            end
        end
        pos = Int64(findmin(assigned_len)[2])
        if total_assigned >= binary_number
            break
        end
    end
    
    ord = sortperm(assigned_len)
    #assign remaining binary variables
    if total_assigned < binary_number
        for i in 1:binary_number
            if has_assigned[i] < 1
                for j in 1 : (threadnum-1)
                    #println(j)
                    if assigned_len[ord[j]] <= assigned_len[ord[j + 1]]
                        has_assigned[i] = j
                        push!(assigned_list[j], i)
                        total_assigned += 1
                        assigned_len[j] += 1
                        break
                    end
                end
                if has_assigned[i] < 1
                    has_assigned[i] = ord[threadnum]
                    push!(assigned_list[ord[threadnum]], i)
                    total_assigned += 1
                    assigned_len[ord[threadnum]] += 1
                end
            end
            #earlier break
            if total_assigned >= binary_number
                break
            end
        end
    end
    
    return assigned_list, has_assigned
end