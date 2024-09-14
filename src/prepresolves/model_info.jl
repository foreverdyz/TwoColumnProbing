using JuMP #use it to read the model
using SparseArrays #constraints info are cached in sparse form


"""
    model_info(filename::String)
This function read the model from the file (named filename), and extract some necessary info for our sequel presolve
"""
function model_info(filename::String)
    
    #read_from_file() is from JuMP
    model = read_from_file(filename)
    #lp_matrix_data() is a new functiom from JuMP (after version 1.7)
    data = lp_matrix_data(model);
 
    """
    suppose the mip is:
        min cTx + c_0
        s.t.  con_lb <= Ax <= con_ub
            var_lb <= x <= var_ub
            x_i in B, i in binary_set
            x_i in Z, i in integer_set
    """
    #cache A to con_matrix
    con_matrix = data.A;
    #cache bounds
    con_ub = data.b_upper;
    con_lb = data.b_lower;
    var_ub = data.x_upper;
    var_lb = data.x_lower;
    #cache a vector for variable name, and we can find the name of the i-th variables
    var_name = data.variables
    
    #cache variables types, and the set includes the orders of these varaibles 
    integer_set = data.integers
    binary_set = data.binaries
    #0 implies continuous, 1 implies integral, and 2 implies binary variable
    var_type = spzeros(length(var_name))
    var_type[integer_set] .= 1
    var_type[binary_set] .= 2
    """
    re-indentify variable types and bounds
    if x is integer and 0<=x<=1, x is binary
    if x is binary, 0<=x<=1 (sometimes its bounds is [-infty, +infty] but it is marked as binary)
    """   
    for index in 1:length(var_ub)
        if var_type[index] > 0
            #check whether an integer variable is a binary variable
            if var_lb[index] == 0 && var_ub[index] == 1
                var_type[index] = 2
            end
            #fix variable bounds of a binary variable
            if var_type[index] > 1
                var_lb[index] = 0
                var_ub[index] = 1
            end
        end
    end
    
    #obj coefficients (c in mip)
    obj_coef = data.c
    #constant in the obj (c_0 in mip)
    obj_constant = data.c_offset
    #obj sense
    if data.sense == MIN_SENSE
        is_min = true
    else
        is_min = false
    end
    
    #then let transfer the way to cache A for later use
    #con_set includes variables' indeces
    con_set = Vector{Int64}[]
    #con_coef includes variables' coefficients
    con_coef = []
    """
    Access terms from the sparse A
    Here we want to get the transpose A, i.e. A^T
    Note that, A' or transpose(A) in Julia does not "really" tranpose A,
    which means the storage structure of A does not change.
    However, we want change the way to cache A so that we can access
    some info fast and with low memory cost.
    """
    I, J, K = findnz(con_matrix);
    push!(I, size(con_matrix)[1])
    push!(J, size(con_matrix)[2])
    push!(K, 0)
    #transpose A
    con_trans = sparse(J, I, K);
    dropzeros!(con_trans);
    #then split variables and coefficients for all costraints
    for i in 1:size(con_matrix)[1]
        local ind, coef = findnz(con_trans[:, i])
        push!(con_set, ind)
        push!(con_coef, coef)
    end
    
    return con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_name, var_type, obj_coef, obj_constant, is_min
end