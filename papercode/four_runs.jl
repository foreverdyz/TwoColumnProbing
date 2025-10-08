#four_runs.jl

using JuMP
using SCIP
using DelimitedFiles
using BangBang

list_1 = open("all_cases.txt") do f
    readlines(f)
end

path1 = "scip_res/4runtime_seed0";

st_term = parse(Int, ARGS[1]);
en_term = parse(Int, ARGS[2]);
memory_GB = 16;
solver_time_limit = 1800;
rand_seed = 0;

for i in st_term:en_term
    filename = list_1[i];
    A = [];
    #Original Model
    model = read_from_file("mipdata/"*filename);
    set_optimizer(model, SCIP.Optimizer);
    set_attribute(model, "parallel/maxnthreads", 1);
    set_attribute(model, "limits/memory", 1024*memory_GB);
    #use the following command for some tests disabling the default probing
    #set_attribute(model, "propagating/probing/maxuseless", false);
    set_attribute(model, "randomization/randomseedshift", rand_seed);
    set_time_limit_sec(model, solver_time_limit);
    optimize!(model);
    if termination_status(model) == TIME_LIMIT
        if primal_status(model) == FEASIBLE_POINT
            rd = [filename[1:end-7], solver_time_limit, MOI.get(model, MOI.ObjectiveBound()), MOI.get(model, MOI.ObjectiveValue())];
        else
            rd = [filename[1:end-7], solver_time_limit, Inf, Inf];
        end
    else
        rd = [filename[1:end-7], solve_time(model), -1, -1];
    end

    #Original Model - probing
    model = read_from_file("mipdata/"*filename);
    set_optimizer(model, SCIP.Optimizer);
    set_attribute(model, "parallel/maxnthreads", 1);
    set_attribute(model, "limits/memory", 1024*memory_GB);
    #use the following command for some tests disabling the default probing
    set_attribute(model, "propagating/probing/maxuseless", 0);
    set_attribute(model, "randomization/randomseedshift", rand_seed);
    set_time_limit_sec(model, solver_time_limit);
    optimize!(model);
    if termination_status(model) == TIME_LIMIT
        if primal_status(model) == FEASIBLE_POINT
            append!!(rd, [solver_time_limit, MOI.get(model, MOI.ObjectiveBound()), MOI.get(model, MOI.ObjectiveValue())]);
        else
            append!!(rd, [solver_time_limit, Inf, Inf]);
        end
    else
        append!!(rd, [solve_time(model), -1, -1]);
    end

    #Reduced Model
    presolve_info = readdlm("presolved_res/presolve_res_16.csv", ',');
    if presolve_info[i, 2]>0
    presolve_time = presolve_info[i, 2];
    model_reduced = read_from_file("presolved_res/presolved_data_16/"*filename);
    set_optimizer(model_reduced, SCIP.Optimizer);
    #use the following command for some tests disabling the default probing
    #set_attribute(model_reduced, "propagating/probing/maxuseless", false);
    set_attribute(model_reduced, "parallel/maxnthreads", 1);
    set_attribute(model_reduced, "limits/memory", 1024*memory_GB);
    set_attribute(model, "randomization/randomseedshift", rand_seed);
    set_time_limit_sec(model_reduced, solver_time_limit - presolve_time);
    optimize!(model_reduced);
    if termination_status(model_reduced) == TIME_LIMIT
        if primal_status(model_reduced) == FEASIBLE_POINT
            append!!(rd, [solver_time_limit, MOI.get(model_reduced, MOI.ObjectiveBound()), MOI.get(model_reduced, MOI.ObjectiveValue())]);
        else
            append!!(rd, [solver_time_limit - presolve_time, Inf, Inf]);
        end
    else
        append!!(rd, [solve_time(model_reduced), -1, -1]);
    end
    else
    append!!(rd, [rd[end-2], rd[end - 1], rd[end]]);
    end

    #Reduced Model - probing
    presolve_info = readdlm("presolved_res/presolve_res_16.csv", ',');
    if presolve_info[i, 2]>0
    presolve_time = presolve_info[i, 2];
    model_reduced = read_from_file("presolved_res/presolved_data_16/"*filename);
    set_optimizer(model_reduced, SCIP.Optimizer);
    #use the following command for some tests disabling the default probing
    set_attribute(model_reduced, "propagating/probing/maxuseless", 0);
    set_attribute(model_reduced, "parallel/maxnthreads", 1);
    set_attribute(model_reduced, "limits/memory", 1024*memory_GB);
    set_attribute(model, "randomization/randomseedshift", rand_seed);
    set_time_limit_sec(model_reduced, solver_time_limit - presolve_time);
    optimize!(model_reduced);
    if termination_status(model_reduced) == TIME_LIMIT
        if primal_status(model_reduced) == FEASIBLE_POINT
            append!!(rd, [solver_time_limit, MOI.get(model_reduced, MOI.ObjectiveBound()), MOI.get(model_reduced, MOI.ObjectiveValue())]);
        else
            append!!(rd, [solver_time_limit - presolve_time, Inf, Inf]);
        end
    else
        append!!(rd, [solve_time(model_reduced), -1, -1]);
    end
    else
    append!!(rd, [rd[end-2], rd[end - 1], rd[end]]);
    end

    open(path1*".csv", "a") do file
        # Write the data
        writedlm(file, [rd], ',')
    end

end

