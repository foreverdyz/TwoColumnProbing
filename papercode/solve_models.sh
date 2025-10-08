#!/bin/bash

for ((i=1; i<=192; i++))
do
 julia --threads 16 solve_models.jl $i $i
done

