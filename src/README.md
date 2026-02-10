# Serial and Parallel Tow-Column Probing for Mixed-Integer Programming

## Presolving Stage

You can presolve an instance with two-column probing method as
```julia
> is_success, time_used, model_reduce, pair_number = presolve(
        filename, size_limit, work_limit, cand_number, max_probe_number, time_limit, threadnum, path
    )
> (is_success > 0) && (write_to_file(model_reduce, "reduced_"*filename*".mps.gz"))
```

is_success = 1 implies two-column probing preprocesses the model; otherwise is_success = 0 means two-column probing has no impact on the model due to time limit or model restriction.

time_uesd is the runtime of two-column probing. 

model_reduce is the presolved model.

pair_number is the number of pairs of binary variables that are probed.

filename is the name of the instance (without .mps or .lp)

size_limit, work_limit, cand_number, max_probe_number, time_limit are parameters introduced in the paper. We set them as 1_000, 200_000_000, 5_000_000, 1_000, 30 in the paper.

threadnum is the number of threads you want to use to parallel (or sequentially) presolve the model. Note that, based on our experiments, we suggest setting it to 1 (serial) or a number >= 8 (parallel).

path is the path to read the file of instance. We read it as path/filename.mps.gz. If your file is not in .mps.gz, you can modify it in line 45 of "presolve.jl"

## Solver Solving Stage

Then you can solve the model with any solver as:
```julia
> model =  read_from_file("reduced_"*filename*".mps.gz")
> set_optimizer(model, solver.optimizer())
#set time limit for the solver, and other features
> optimize!(model)
```

## Combine presolving and SCIP Solving

You can also presolve the model and solve the reduced model with SCIP, like our paper:

```julia
> solver_time, presolve_time, dual_bound, primal_bound =  solve_model(
        filename, size_limit, work_limit, cand_number, max_probe_number, 
        time_limit, threadnum, path, solve_org, disable_probing,
        solver_threads, memory_GB, solver_time_limit
    )
```

solve_org = 1 implies you are solving the original model without two-column probing; otherwise, solve_org = 0 means you are presolving with two-column probing and then solve it with SCIP.

disable_probing = 1 implies disabling the classical probing method in SCIP; otherwise, disable_probing = 0 to keep the original probing.

solver_threads is the number of threads used in SCIP, which in the current version (SCIP v0.11.14), it only impacts the presolving method.

memory_GB is the limit of memory used by the solver, which we suggest setting to your PC's memory.

solver_timie_limit is the time limit of the solver, which is set as 3600 in default.





