#classify_constraints.jl

using SparseArrays

"""
    Here we want to find:
    first, set packing constraint;
    second, other constraints with binary variables;
    I includes indeces for all mix-binary constraint, J inlcudes all clque indeces; B inlcudes binary term for I, and S for J
"""
function classify_constraints(con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64})
    #there are m constraints
    m = length(con_set)
    """
        Initialize three index sets: I for pure binary constraints, and J for mixed binary constraints
    """
    I = Int64[]
    J = Int64[]
    S = Vector{Int64}[]
    B = Vector{Int64}[]
    for i in 1:m
        #determine whether this constraint is pure binary, mixed binary, or non-binary
        #id = 0 (non-binary), 1 (others), 2 (potential set packing)
        id, bin_set = _con_class(con_set[i], con_coef[i], var_type, org_to_bin)
        if id == 1
           push!(I, i)
           push!(B, bin_set)
        elseif id == 2
            local id_first = false
            local id_second = false
            #treat one constraints as two, >=  and <=
            if con_ub[i] != Inf
                s, id_first = _check_set_pack(con_set[i], con_coef[i], con_ub[i], org_to_bin)
                if id_first
                    push!(S, s)
                    push!(J, i)
                else
                    push!(I, i)
                    push!(B, bin_set)
                end
            end
            if con_lb[i] != -Inf
                #if the constraint has been identified as set packing constraint
                #just denote the other side of the constraint as a mix-binary constriant.
                if id_first
                    push!(I, i)
                    push!(B, bin_set)
                else
                    s, id_second = _check_set_pack(con_set[i], -con_coef[i], -con_lb[i], org_to_bin)
                    if id_second
                        push!(S, s)
                        push!(J, i) 
                    end
                end
            end
        end
    end
    return I, S, B, J
end

"""
    con_class()
    this classification function determines whether one constraint is pure binary, mixed binary, or non-binary
"""
function _con_class(con::Vector{Int64}, coef::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64})
    """
        return 0 (non-binary), 1 (others), 2 (potential set packing) 
    """
    let
        #s cache all binary variables
        local s = Int64[]
        #is_set implies whether this constraint is a set packing constraint
        local is_set = true
        for (i, j) in enumerate(con)
            (var_type[j] > 1) && (push!(s, org_to_bin[j]))
            if is_set
                (abs(coef[i]) != 1) && (is_set = false)
            end
        end
        if is_set
            #we want all variables are binary variables
            (length(s) == length(con)) && (return 2, s)
        end
        #no binary variables
        (length(s) == 0) && (return 0, s)
        return 1, s
    end
end

function _check_set_pack(con::Vector{Int64}, coef::AbstractVector, b::Real, org_to_bin::Dict{Int64, Int64})
    binary_length = length(org_to_bin)
    #build an empty expended constraint
    s = zeros(length(con))
    is_set_pack = true
    for (i, j) in enumerate(con)
        if coef[i] > 0
            s[i] = org_to_bin[j]
        else
            s[i] = org_to_bin[j] + binary_length
            b += 1
        end
    end
    (b != 1) && (is_set_pack = false)
    return s, is_set_pack
end