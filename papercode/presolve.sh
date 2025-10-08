#!/bin/bash

julia --threads 1 presolve_models.jl 1 1 192
julia --threads 2 presolve_models.jl 2 1 192
julia --threads 4 presolve_models.jl 4 1 192
julia --threads 8 presolve_models.jl 8 1 192
julia --threads 16 presolve_models.jl 16 1 192
