# Serial and Parallel Tow-Column Probing for Mixed-Integer Programming

## Presolving Stage

You can presolve an instance with two-column probing method as
```julia
> is_success, time_used, model_reduce, pair_number = presolve(filename, size_limit, work_limit, cand_number, max_probe_number, time_limit, threadnum, path)
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

## SCIP Solving Stage


## Disable SCIP's Probing




