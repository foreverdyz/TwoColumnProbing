# TwoColumnProbing

This archive is under the [GNU General Public License v3.0](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in [Serial and Parallel Two-Column Probing for Mixed-Integer Programming](https://arxiv.org/abs/2408.16927) by Yongzheng Dai and Chen Chen.

## Something Important!!!

Please install an alpha version of JuMP to reproduce our result! The prototype can be obtained by:

```julia
pkg> add JuMP  #Press ']' to enter the Pkg REPL mode.

pkg> add MathOptInterface#master
``` 

## Language

The code is written in the [Julia programming language](https://julialang.org). Please visit the website for an installation guide. 

## Required packages 

To run the content of scripts in [src](/src), the following Julia packages are required: SCIP, BangBang, Random, SparseArrays, JuMP.

Yon can check whether you have installed these required packages by running [package_check.jl (in src)](/src/packge_check.jl)
```julia
julia> include("package_check.jl")
```

To install a package, simply run

```julia
pkg> add PACKAGE_NAME    # Press ']' to enter the Pkg REPL mode.
```

## Tutorial

See [README.md (in src)](/src/README.md) for tutorial and examples.

