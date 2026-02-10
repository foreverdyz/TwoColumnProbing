# Serial and Parallel Tow-Column Probing for Mixed-Integer Programming

## Test Cases
In our manuscript, there are 190 cases and some removed cases because of memory exceeding. Now we solved the memroy issue and there are 192 instances saved
in "all_cases.txt". However, when using our two-column probing, there still be some cases are not impacted by our method, e.g. all pairs of binary variables
have had conflicts between each other. As a result, in "presolved\_res/presolve\_res\_threadnum.csv", some instances are detected 0 pair of variables. In this
case, we will not test their runtime with SCIP.

## Software Versions

Julia 1.8.5 

JuMP 1.28 

SCIP v0.11.14

## Presolving Stage (Experiments for Sec. 5.3)

Run "Nohup bash presolve.sh" to generate presolved results, which will be saved in Folder "presolved\_res". The presolving runtime results correspond to Sec. 5.3: Table 1 and Figure 2.

Note that: sometimes the program maybe killed. The potential reason might be continuously running too many instances in one program and the server kills it.
If it happens, you can check in which case the program is killed in "presolved\_res/presolve\_res\_threadnum.csv" where threadnum is the number of threads.
For example, if it is killed in threadnum = 8 and the 100th instance, you can open "presolve.sh" and change the file as 

```julia
#comment threadnum = 1,2,4
julia --threads 8 presolve_models.sh 8 100 192
julia --threads 16 presolve_models.sh 16 1 192
```


## SCIP Solving Stage (Experiments for Sec. 5.4)

Run "Nohup bash solve_models.sh" to solve all reduced models with SCIP and all results will be saved in Folder "SCIP\_res". The SCIP runtime results correspond to Sec. 5.4: Table 2 and Table 3.

Note that, due to the same reason, the program is sometimes killed. You can use the similar method to solve this question by checking "SCIP\_res/runtime\_1\_192.csv".

## Disable SCIP's Probing (Experiments for Sec. 5.5)

There are some experiments in our manuscript that diasbled SCIP's original probing method. You can do so by checking script "solve_models.jl". There are some lines
we left comments for how to disable SCIP's probing. (Just need to uncomment some lines). The SCIP runtime results correspond to Sec. 5.5: Table 4 and Table 5.

## Results Comp (How we calculate all 1-shifted geometric mean results.)

Runtime results will be saved in "SCIP\_res/runtime\_1\_192.csv".

```julia
using StatsBase
#suppose list a caches runtime results for org scip, b caches runtimes of one presolved method, and c is corresponding presolving time
I, J = Int64[], Int64[]
#I includes instances that at least one method can solve it within timelimit
#J includes instances taht at least one method cannot solve it within timelimit
for i in 1:length(a)
    if a[i] >= 3600
        push!(J, i)
	if b[i] + c[i] >= 3600
	    push!(I, i)
    	end
    else
        if b[i] + c[i] >= 3600
           push!(J, i)
        end
    end
end

runtime1, runtime2 = geomean(a[I].+10), geomean(b[I].+10)

#suppose list d, e, f, g caches upper and lower bounds for org scip and one presolved method
gap1 = geomean([abs(d[i] - e[i])/max(abs[d[i]], abs[e[i]]) for i in J].+1)
gap2 = geomean([abs(f[i] - g[i])/max(abs[f[i]], abs[g[i]]) for i in J].+1)
```

# Experiments in Sec. 5.6

## Compare 4 methods with different random seeds

Run "Nohup bash four_runs.sh" to solve all models with 4 methods, default SCIP, default SCIP disable Probing, Two-Column Probing + default SCIP, Two-Column Probing + default SCIP disable Probing. All results will be saved in Folder "SCIP\_res". The SCIP runtime results correspond to Sec. 5.6, Table 7.

If you want to change the random seed, please open file "four_runs.jl", and change line 12
```julia
path1 = "scip_res/4runtime_seed0"
```
and line 18
```julia
rand_seed = 0
```
to any other seeds, for example path1 = "scip_res/4runtime_seed1" and rand_seed = 1. We use 0, 1, 2, 3, 4.



