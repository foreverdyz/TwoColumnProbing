#solve_model.jl

include("presolve.jl");

using SCIP

"""
    solve_model()
Preosolve model with two-column probing and solve it with scip.
Intput:
filename: the name of .mps file
size_limit: the maximum size of constraints we considered in building coupling matrix
work_limit: the maximum number of elements processed in building coupling matrix
cand_number: the number of candidates we considered in probing order
max_probe_number: the maximum number of pairs of variables that we may probe (the upper bound of number of probing iterations)
time_limit: the time limit for our presolving
threadnum: the number of threads we used
path: the dir saving the original models
solve_org: > 0 implies solving the original model
disable_probing: >0 will disabling probing in SCIP
solver_threads: number of threads used by the solver
memory_GB: number of GBs of memory used by the solver
solver_time_limit: solver's time limit. Note that, we will remove the presolving time. 
"""
function solve_model(
        filename::String, size_limit::Int64, work_limit::Int64, 
        cand_number::Int64, max_probe_number::Int64, 
        time_limit::Real, threadnum::Int64, path::String, solve_org::Int64,
        disable_probing::Int64, solver_threads::Int64, memory_GB::Real, solver_time_limit::Real = 3600
    )
    if solve_org > 0
        model = read_from_file(path*"/"*filename*".mps");
        set_optimizer(model, SCIP.Optimizer);
        set_attribute(model, "parallel/maxnthreads", solver_threads);
        if disable_probing > 0
            set_attribute(model, "propagating/probing/maxuseless", false);
        end
        set_attribute(model, "limits/memory", 1024*memory_GB)
        set_time_limit_sec(model, solver_time_limit);
        set_silent(model);
        optimize!(model);
        if termination_status(model) == TIME_LIMIT
            @warn "Solver do not solve the problem within time_limit"
            if primal_status(model) == FEASIBLE_POINT
                @info "Solver found an feasible point, will return solver_time_limit, 0, primal bound, dual bound"
                return solver_time_limit, 0, MOI.get(model, MOI.ObjectiveBound()), MOI.get(model, MOI.ObjectiveValue())
            else
                @warn "Solver did not find any feasible point, will return solver_time_limit, 0, Inf, Inf"
                return solver_time_limit, 0, Inf, Inf
            end
        else
            @info "Solver solved the model. Will return solver time, 0, -1, -1"
            return solve_time(model), 0, -1, -1
        end
    else
        id, time_used, model_reduce, _ =  presolve(
            filename, size_limit, work_limit, cand_number, max_probe_number, time_limit, threadnum, path
        );
        if id < 0
            @warn "Runtime Exceeded, please set a longer runtime. Will retun (-1, -1, -1, -1)."
            return -1, -1, -1, -1
        elseif id > 0
            set_optimizer(model_reduce, SCIP.Optimizer);
            set_attribute(model_reduce, "parallel/maxnthreads", solver_threads);
            if disable_probing > 0
                set_attribute(model_reduce, "propagating/probing/maxuseless", false);
            end
            set_attribute(model_reduce, "limits/memory", 1024*memory_GB)
            set_time_limit_sec(model_reduce, solver_time_limit - time_used);
            set_silent(model_reduce);
            optimize!(model_reduce);
            if termination_status(model_reduce) == TIME_LIMIT
                @warn "Solver do not solve the problem within time_limit"
                if primal_status(model_reduce) == FEASIBLE_POINT
                    @info "Solver found an feasible point, will return solver_time_limit, 0, primal bound, dual bound"
                    return solver_time_limit, 0, MOI.get(model_reduce, MOI.ObjectiveBound()), MOI.get(model_reduce, MOI.ObjectiveValue())
                else
                    @warn "Solver did not find any feasible point, will return solver_time_limit, 0, Inf, Inf"
                    return solver_time_limit, 0, Inf, Inf
                end
            else
                @info "Solver solved the model. Will return solver time, presolving time, -1, -1"
                return solve_time(model_reduce), time_used, -1, -1
            end
        else
            @warn "Two-Column probing has no impact to the model. Will return (0, 0, -1, -1)."
            return 0, 0, -1, -1
        end
    end
end