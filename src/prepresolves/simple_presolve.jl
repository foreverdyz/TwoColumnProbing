#simple_presolve.jl

using SparseArrays

"""
The simple presolve includes one round of singleton removing, one-row variables bounds strengthenin, and binary re-identification
There are following reasons to do so:
first, the probing method is very sensitive to the variables' bounds, and some initial bounds of variables are too loose;
second, these three methods are very easy to implement and do not have a big impact to the runtime (see the gurobi ijoc paper)
"""
function simple_presolve(A::AbstractSparseArray, con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #identify empties and singletons and update variables' bounds
    I, var_ub, var_lb = _singleton_remove(con_set, con_coef, con_ub, con_lb, var_ub, var_lb)
    con_ub, con_lb, var_ub, var_lb = _bounds_strengthen(con_set, con_coef, con_ub, con_lb, var_ub, var_lb, I, var_type)
    #remove empties and singletons
    A = A[I, :]
    con_set = con_set[I]
    con_coef = con_coef[I]
    #re-indentify binary variables and get maps between original variables and binary variables
    var_type, org_to_bin, bin_to_org = _iden_bin_list(var_ub, var_lb, var_type)
    return A, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, org_to_bin, bin_to_org
end

"""
    singleton_remove()
    this function removes singleton from the constraints matrix and update variable bounds 
"""
function _singleton_remove(con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector)
    #get m constraints
    m = length(con_set)
    #index set for non-singleton constraint
    I = Int64[]
    #check all constraints
    for j in 1:m
        local con = con_set[j]
        if length(con) == 1 #which implies len(c) is 1
            #i is the index of the only variable in the j-th constraints
            local i = con[1]
            #get the coefficient
            local c = con_coef[j][1]
            #then we can strengthen variables' bounds
            if c > 0
                (con_ub[j] != Inf) && (var_ub[i] = min(var_ub[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_lb[i] = max(var_lb[i], con_lb[j]/c))
            else
                (con_ub[j] != Inf) && (var_lb[i] = max(var_lb[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_ub[i] = min(var_ub[i], con_lb[j]/c))
            end
        elseif length(con) > 1 #remove empty constraints
            #j is not a singleton constraint; thus push it to the index set I
            #we also remove too longer constriants
            (length(con) <= 1_000_000) && (push!(I, j))
        end
    end
    #output the index set I and updated varaibles' bounds
    return I, var_ub, var_lb
end

"""
    bounds_stren()
    this function propagates a simple variable bounds strengthening
"""
function _bounds_strengthen(con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, I::Vector{Int64}, var_type::AbstractVector)
    for j in I
        """
            the logic is as follows:
            1, <= and >= will be treated as one inequality
            2, == will be treated as two inequalities separately
        """
        (con_ub[j] != Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(con_set[j], con_coef[j], con_ub[j], var_ub, var_lb, var_type))
        (con_lb[j] != -Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(con_set[j], -con_coef[j], -con_lb[j], var_ub, var_lb, var_type))
    end
    #output constraints matrix and bounds with only constraints from I and updated variables bounds
    return con_ub[I], con_lb[I], var_ub, var_lb
end

function _unique_ub_con_bound_strengthen(con::Vector{Int64}, coef::AbstractVector, b::Real,var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #ub means we want lower bouond for every variable; we use delta_list to record each lower bound in this constraint
    delta_list = zeros(length(con))
    for (i, j) in enumerate(con)
        #calculate lower bound for each variables
        (coef[i] > 0) ? (delta_list[i] = coef[i]*var_lb[j]) : (delta_list[i] = coef[i]*var_ub[j])
    end
    delta = sum(delta_list)
    (delta == -Inf) && (return var_ub, var_lb)
    for (i, j) in enumerate(con)
        #the new potential bound
        local x = (b - delta + delta_list[i])/coef[i]
        #====
            Here we do not need to udpate the delta, since the updated variable bound is not to do with the delta!
        ====#
        if coef[i] > 0
            if var_ub[j] > x + 0.00001
                (var_type[j] > 0) ? (var_ub[j] = floor(Int64, x)) : (var_ub[j] = x)
            end
        else
            if var_lb[j] < x - 0.00001
                (var_type[j] > 0) ? (var_lb[j] = ceil(Int64, x)) : (var_lb[j] = x)
            end
        end
    end
    #output
    return var_ub, var_lb
end

function _unique_lb_con_bound_strengthen(con::Vector{Int64}, coef::AbstractVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    delta_list = zeros(length(con))
    for (i, j) in enumerate(con)
        #calculate lower bound for each variables
        (coef[i] > 0) ? (delta_list[i] = coef[i]*var_ub[j]) : (delta_list[i] = coef[i]*var_lb[j])
    end
    delta = sum(delta_list)
    (delta == Inf) && (return var_ub, var_lb)
    for (i, j) in enumerate(con)
        #the new potential variable bound
        local x = (b - delta + delta_list[i])/coef[i]
         #====
            Here we do not need to udpate the delta, since the updated variable bound is not to do with the delta!
        ====#
        if coef[i] < 0
            if var_ub[j] > x + 0.00001
                (var_type[j] > 0) ? (var_ub[j] = floor(Int64, x)) : (var_ub[j] = x)
            end
        else
            if var_lb[j] < x - 0.00001
                (var_type[j] > 0) ? (var_lb[j] = ceil(Int64, x)) : (var_lb[j] = x)
            end
        end
    end
    #output
    return var_ub, var_lb
end

"""
    iden_bin_list()
    this function implements two operations: first, re-indentify (potential) binary variables from integer variables; second,
build a binary map (including bin to org and org to bin) to help us find a binary variable from original variables and find
the original variable for a binary variable.
"""
function _iden_bin_list(var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    org_to_bin = Dict{Int64, Int64}()
    bin_to_org = Dict{Int64, Int64}()
    index = 0
    for (i, v) in enumerate(var_type)
        if v == 1
            if var_ub[i] == 1 && var_lb[i] == 0
                var_type[i] = 2
                index += 1
                org_to_bin[i] = index
                bin_to_org[index] = i
            end
        elseif v == 2
            index += 1
            org_to_bin[i] = index
            bin_to_org[index] = i
        end
    end
    return var_type, org_to_bin, bin_to_org
end