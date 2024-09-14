#get_presolved_model.jl

include("presolve.jl");

"""
    get_presolved_model()
Will presolve the model with two-column probing method, then write the model to a .mps file.
Intput:
filename: the name of .mps file
size_limit: the maximum size of constraints we considered in building coupling matrix
work_limit: the maximum number of elements processed in building coupling matrix
cand_number: the number of candidates we considered in probing order
max_probe_number: the maximum number of pairs of variables that we may probe (the upper bound of number of probing iterations)
time_limit: the time limit for our presolving
threadnum: the number of threads we used
org_path: the dir saving the original model
save_path: the dir you want to save the file
info_id: > 0  implies we will provide some important info
"""
function get_presolved_model(
        filename::String, size_limit::Int64, work_limit::Int64, 
        cand_number::Int64, max_probe_number::Int64, 
        time_limit::Real, threadnum::Int64, org_path::String, save_path::String,
        info_id::Int64 = 1
    )
    info0 = "Writing the model to a file may take extra runtime.\n"
    info1 = "We recommend to solve it in Julia directly without writing down.\n"
    info2 = "Julia JuMP may read .mps file and preprocessing it a bit, which is different to other solvers or languages.\n"
    info3 = "Thus, comparing org and reduced models in other env may return unfair results.\n"
    info4 = "Please set solvers in JuMP to org and redcued models and solve them.\n"
    info5 = "You can turn off these info by setting info_id = 0"
    (info_id > 0) && (@info info0*info1*info2*info3*info4*info5)
    id, time_used, model_reduce, _ = presolve(
        filename, size_limit, work_limit, cand_number, max_probe_number, time_limit, threadnum, org_path
    );
    if id < 0
        @warn "Runtime Exceeded, please set a longer runtime. Will retun -1."
        return -1
    elseif id > 0
        write_to_file(model_reduce, save_path*"/"*filename*"_"*string(threadnum)*"_reduced.mps");
        @info "Save reduced model to "*save_path*"/"*filename*"_"*string(threadnum)*"_reduced.mps. Will return presolving time."
        return time_used
    else
        @warn "Two-Column probing has no impact to the model. Will return 0."
        return 0
    end
end