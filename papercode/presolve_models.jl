#presolve_models.jl

include("presolve.jl");

using JuMP
using DelimitedFiles

list_1 = open("all_cases.txt") do f
    readlines(f)
end;

threadnum = parse(Int, ARGS[1]);
st_term = parse(Int, ARGS[2]);
en_term = parse(Int, ARGS[3]);

A = [];

#precompile presolve()
presolve(
    list_1[1][1:end-7], 1_000, 200_000_000, 5_000_000,
    1_000, 30, threadnum, "mipdata"
);

for i in st_term:en_term
    filename = list_1[i][1:end-7];
    print(i, " ", filename);
    id, time_used, model, pair_number = presolve(
        filename, 1_000, 200_000_000, 5_000_000, 
        1_000, 30, threadnum, "mipdata" 
    );
    if id > 0
        write_to_file(model, "presolved_res/presolved_data_"*string(threadnum)*"/"*filename*".mps.gz");
	open("presolved_res/presolve_res_"*string(threadnum)*".csv", "a") do file
            # Write the data
            writedlm(file, [[filename, time_used, pair_number]], ',')
        end
    else
	open("presolved_res/presolve_res_"*string(threadnum)*".csv", "a") do file
            # Write the data
            writedlm(file, [[filename, 0]], ',')
        end
    end
    println(" Done");
    model = nothing;
    GC.gc();
end

