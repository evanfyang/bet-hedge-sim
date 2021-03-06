#!/usr/bin/env julia

## vc_submit_grid.jl
##
## Script for submiting params to run 
## simulations in batch via SLURM;
## Adapted from vc_submit_grid.jl and parallel_simulate_population.jl by Taylor Kessinger

using JLD
using ArgParse
using Dates
using Distributed
using DataFrames
using CSV
import JSON
include("./read_parameters.jl")


function main(args)
    s = ArgParseSettings(description = "Run bet hedge simulation across multiple cores")
        @add_arg_table s begin
            "--ncpus"
                arg_type=Int64
                default=max(round(Int,Sys.CPU_THREADS/2), 1)
            "--input"
                default=nothing
            "--output"
                arg_type=String
                default="test"
        end
        parsed_args = parsed_args(args, s)

    jobids = []
    cmd_strings = []

    date = Dates.format(Dates.Date(Dates.now()), "yyyymmdd")

    defpars = Dict{String,Any}([
        "num_gen"           => Dict("value" => 10^5, "type" => Int64),
        "pop_size"           => Dict("value" => 1e-5, "type" => Int64),
        "env_switch_rate"   => Dict("value" => 1e-3, "type" => Int64),
        "init_switch_rate"  => Dict("value" => 1e-2, "type" => Float64),
        "mutation_rate"     => Dict("value" => 1.0,  "type" => Float64),
        "s1"                => Dict("value" => 0.01, "type" => Float64),
        "s2"                => Dict("value" => 1e-5, "type" => Float64),
    ])
    
    pars = read_parameters(defpars, parsed_args["input"])
    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    runs_per_param_comb = 10

    basename = date*"_grid_tunnel"

    for p in 1:nsets
        for runno in 1:runs_per_param_comb

            simstr = "../src/simulation.jl"

            # create filname base from basename and run number
            numstr = lpad(string(p), length(string(nsets)), "0")
            runstr = lpad(string(runno), 2, "0")
            filebase = basename * "_" * numstr * "_" * runstr

            # add parameter values as command line options
            simstr *= " "
            simstr *= join(map( (x) -> "--" * x[1] * " " * string(x[2]), zip(pars, parsets[p])), " ")
            simstr *= " --file " * filebase * ".jld"

            # run the scripts via sbatch
            # the "split" wizardry is needed due to the fine distinction between Cmd and String types
            println(`sbatch $(split(simstr))`)
            jobout = String(read(`sbatch $(split(simstr))`))

            #jobout = String(read(`sbatch $simstr`))

            push!(cmd_strings, `sbatch $(split(simstr))`)
            push!(jobids, match(r"[0-9]+", jobout).match)
        end
    end

    jid_file = "$(date)_grid_sims_jobids.dat"
    cmd_file = "$(date)_grid_sims_commands.dat"

    writedlm(jid_file, [[(jobids[p], p, runno) for runno in runs_per_param_comb] for p in 1:nsets])
    println("\njob IDs written to file $jid_file\n")
    writedlm(cmd_file, [(jobids[x], cmd_strings[x]) for x in length(jobids)])
    println("submitted commands written to file $cmd_file\n")

main(ARGS)