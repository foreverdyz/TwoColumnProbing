#!/bin/bash

for ((i=1; i<=192; i++))
do
 julia --threads 16 four_runs.jl $i $i
done

