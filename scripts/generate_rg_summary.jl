#!/usr/bin/env julia
# =============================================================================
# generate_rg_summary.jl - Generate Rg summary in block-rg-all.csv format
# =============================================================================
#
# Usage:
#   julia --project=. scripts/generate_rg_summary.jl \
#       --input-dir _data/Block-N200-passive2M-active3M-perrun/csv \
#       --output _data/Block-N200-passive2M-active3M-perrun/rg_summary.csv \
#       --tail-percent 20
#
# Output format: N,x,y,y_err
#   - N: number of monomers
#   - x: n_active/N (fraction of active monomers)
#   - y: Rg/Rg(x=0) (normalized by passive reference)
#   - y_err: standard error of y
#
# Note: For x=0 reference, prefers active_rg.csv files with data, falls back to
#       thermal_rg.csv when active files are empty (consistent with plot_rg_normalized.jl)
#
# =============================================================================

using ArgParse
using CSV
using DataFrames
using Statistics
using Glob

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate Rg summary in normalized format (N, x, y, y_err)",
        prog = "generate_rg_summary.jl"
    )

    @add_arg_table! s begin
        "--input-dir"
            help = "Directory containing CSV files"
            arg_type = String
            required = true
        "--output"
            help = "Output CSV file path"
            arg_type = String
            required = true
        "--tail-percent"
            help = "Percentage of tail data to use for averaging"
            arg_type = Float64
            default = 20.0
        "--verbose"
            help = "Print detailed progress"
            action = :store_true
    end

    return parse_args(s)
end

"""
Parse sweep filename to extract N and n_active.
Returns (n_monomers, n_active, is_thermal).
"""
function parse_filename(filename::String)
    basename_str = basename(filename)

    # Check if it's a thermal file (passive equilibration)
    is_thermal = occursin("thermal_rg.csv", basename_str)

    # Match pattern: single_N{n}_Nact{nact}_..._rg.csv
    m = match(r"single_N(\d+)_Nact(\d+)_", basename_str)
    if m !== nothing
        return (
            n_monomers = parse(Int, m[1]),
            n_active = parse(Int, m[2]),
            is_thermal = is_thermal
        )
    end

    return nothing
end

"""
Compute mean Rg from the tail portion of a time series.
"""
function compute_tail_mean(filepath::String; tail_percent::Float64=20.0)
    df = CSV.read(filepath, DataFrame)

    n = nrow(df)
    if n == 0
        return NaN
    end

    tail_start = max(1, n - round(Int, n * tail_percent / 100) + 1)
    tail_data = df[tail_start:end, :]

    # Get Rg column
    if "Rg" in names(df)
        rg_values = collect(skipmissing(tail_data.Rg))
        valid = filter(!isnan, rg_values)
        return isempty(valid) ? NaN : mean(valid)
    else
        error("No 'Rg' column found in $filepath")
    end
end

function main()
    args = parse_commandline()

    input_dir = args["input-dir"]
    output_path = args["output"]
    tail_percent = args["tail-percent"]
    verbose = args["verbose"]

    # Find active_rg files (for all x, including x=0 if they have data)
    active_pattern = joinpath(input_dir, "single_N*_Nact*_*_active_rg.csv")
    active_files = glob(active_pattern)

    # Find thermal_rg files (fallback for x=0 when active files are empty)
    thermal_pattern = joinpath(input_dir, "single_N*_Nact0_*_thermal_rg.csv")
    thermal_files = glob(thermal_pattern)

    println("Found $(length(active_files)) active_rg files")
    println("Found $(length(thermal_files)) thermal_rg files (fallback for n_active=0)")

    # Group files by (N, n_active) and compute tail means
    data = Dict{Tuple{Int,Int}, Vector{Float64}}()

    # Process active files first (including n_active=0 if they have data)
    for file in active_files
        params = parse_filename(file)
        if params === nothing
            verbose && println("  Skipping (unparseable): $(basename(file))")
            continue
        end

        key = (params.n_monomers, params.n_active)
        rg_mean = compute_tail_mean(file; tail_percent=tail_percent)

        if isnan(rg_mean)
            verbose && println("  Skipping (no data): $(basename(file))")
            continue
        end

        if !haskey(data, key)
            data[key] = Float64[]
        end
        push!(data[key], rg_mean)

        verbose && println("  $(basename(file)) -> Rg_tail = $rg_mean")
    end

    # Process thermal files as fallback for n_active=0 when no active data exists
    for file in thermal_files
        params = parse_filename(file)
        if params === nothing
            verbose && println("  Skipping (unparseable): $(basename(file))")
            continue
        end

        key = (params.n_monomers, 0)  # Force n_active=0 for thermal

        # Only use thermal as fallback if no active data exists for this (N, 0) key
        if haskey(data, key) && !isempty(data[key])
            verbose && println("  Skipping (active data exists): $(basename(file))")
            continue
        end

        rg_mean = compute_tail_mean(file; tail_percent=tail_percent)

        if isnan(rg_mean)
            verbose && println("  Skipping (no data): $(basename(file))")
            continue
        end

        if !haskey(data, key)
            data[key] = Float64[]
        end
        push!(data[key], rg_mean)

        verbose && println("  $(basename(file)) -> Rg_tail = $rg_mean (thermal fallback)")
    end

    # Compute statistics for each (N, n_active) group
    results = DataFrame(N=Float64[], x=Float64[], y=Float64[], y_err=Float64[])

    # Sort keys by (N, n_active)
    sorted_keys = sort(collect(keys(data)))

    # First, find Rg(x=0) for each N (for normalization)
    rg_passive = Dict{Int, Float64}()
    for (n, nact) in sorted_keys
        if nact == 0
            values = filter(!isnan, data[(n, nact)])
            if !isempty(values)
                rg_passive[n] = mean(values)
                println("  Passive reference for N=$n: Rg0 = $(rg_passive[n]) ($(length(values)) runs)")
            end
        end
    end

    # Now compute normalized results
    for (n, nact) in sorted_keys
        values = filter(!isnan, data[(n, nact)])
        if isempty(values)
            continue
        end

        rg_mean = mean(values)
        rg_sem = length(values) > 1 ? std(values) / sqrt(length(values)) : 0.0

        # Normalize by passive reference
        if !haskey(rg_passive, n)
            @warn "No passive (n_active=0) reference for N=$n, using absolute values"
            rg0 = 1.0
        else
            rg0 = rg_passive[n]
        end

        x = nact / n
        y = rg_mean / rg0
        y_err = rg_sem / rg0  # Error propagation (assuming rg0 has negligible error)

        push!(results, (N=Float64(n), x=x, y=y, y_err=y_err))
    end

    # Sort by N, then by x
    sort!(results, [:N, :x])

    # Save results
    output_dir = dirname(output_path)
    if !isempty(output_dir) && !isdir(output_dir)
        mkpath(output_dir)
    end

    CSV.write(output_path, results)
    println("\nSaved to: $output_path")
    println("\nSummary:")
    println(results)
end

main()
