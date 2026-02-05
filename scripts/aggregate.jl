#!/usr/bin/env julia
# =============================================================================
# aggregate.jl - Aggregate metrics from multiple simulation runs
# =============================================================================
#
# Usage:
#   # Aggregate files matching a pattern
#   julia --project=. scripts/aggregate.jl \
#       --pattern "_data/csv/single_100_50_0.0_5.0_*_active_msd.csv" \
#       --output "_data/csv/aggregated_msd.csv"
#
#   # Aggregate all metrics for a base name
#   julia --project=. scripts/aggregate.jl \
#       --base "single_100_50_0.0_5.0" \
#       --phase active \
#       --input-dir "_data/csv" \
#       --output-dir "_data/csv/aggregated"
#
# =============================================================================

using ArgParse
using CSV
using DataFrames
using Statistics
using Glob

"""
Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings(
        description = "Aggregate metrics from multiple simulation runs",
        prog = "aggregate.jl"
    )

    @add_arg_table! s begin
        "--pattern"
            help = "Glob pattern to match input CSV files"
            arg_type = String
            default = ""
        "--base"
            help = "Base name to find all metrics (alternative to --pattern)"
            arg_type = String
            default = ""
        "--phase"
            help = "Phase to aggregate: active or thermal"
            arg_type = String
            default = "active"
        "--output"
            help = "Output file path (for single pattern mode)"
            arg_type = String
            default = ""
        "--input-dir"
            help = "Input directory for CSV files"
            arg_type = String
            default = "_data/csv"
        "--output-dir"
            help = "Output directory for aggregated files (for --base mode)"
            arg_type = String
            default = "_data/csv/aggregated"
        "--interpolate"
            help = "Interpolate to common time points if times don't match exactly"
            action = :store_true
    end

    return parse_args(s)
end

"""
Find files matching a glob pattern
"""
function find_matching_files(pattern::String)
    files = glob(pattern)
    if isempty(files)
        error("No files found matching pattern: $pattern")
    end
    return sort(files)
end

"""
Find all metric files for a base name and phase
"""
function find_metric_files(base::String, phase::String, input_dir::String)
    metrics = Dict{String, Vector{String}}()

    # Common metric suffixes
    metric_types = ["msd", "rg"]

    for metric in metric_types
        pattern = joinpath(input_dir, "$(base)_*_$(phase)_$(metric).csv")
        files = glob(pattern)
        if !isempty(files)
            metrics[metric] = sort(files)
        end
    end

    if isempty(metrics)
        error("No metric files found for base='$base', phase='$phase' in '$input_dir'")
    end

    return metrics
end

"""
Read and validate CSV files, ensuring they have compatible structure
"""
function read_csv_files(files::Vector{String})
    dfs = DataFrame[]

    for file in files
        df = CSV.read(file, DataFrame)
        push!(dfs, df)
    end

    # Check that all files have the same columns
    ref_cols = names(dfs[1])
    for (i, df) in enumerate(dfs)
        if names(df) != ref_cols
            @warn "File $(files[i]) has different columns than $(files[1])"
        end
    end

    return dfs
end

"""
Get the time column name from a dataframe
"""
function get_time_column(df::DataFrame)
    possible_names = ["lag_time", "time", "t", "step"]
    for name in possible_names
        if name in names(df)
            return name
        end
    end
    # Default to first column
    return names(df)[1]
end

"""
Aggregate dataframes by computing mean, std, and SEM for each column
"""
function aggregate_dataframes(dfs::Vector{DataFrame}; interpolate::Bool=false)
    n_runs = length(dfs)
    time_col = get_time_column(dfs[1])

    # Get all numeric columns except time
    numeric_cols = [name for name in names(dfs[1])
                    if name != time_col && eltype(dfs[1][!, name]) <: Number]

    if interpolate
        # Find common time range and interpolate
        return aggregate_with_interpolation(dfs, time_col, numeric_cols, n_runs)
    else
        # Assume all files have same time points
        return aggregate_exact(dfs, time_col, numeric_cols, n_runs)
    end
end

"""
Aggregate assuming exact time point matching
"""
function aggregate_exact(dfs::Vector{DataFrame}, time_col::String,
                         numeric_cols::Vector{String}, n_runs::Int)
    # Use time points from first file
    times = dfs[1][!, time_col]
    n_points = length(times)

    # Check that all files have same number of points
    for (i, df) in enumerate(dfs)
        if nrow(df) != n_points
            @warn "File $i has $(nrow(df)) rows, expected $n_points. Will use minimum."
        end
    end

    # Use minimum number of points across all files
    n_points = minimum(nrow.(dfs))
    times = times[1:n_points]

    # Create result dataframe
    result = DataFrame()
    result[!, time_col] = times

    for col in numeric_cols
        # Collect values from all runs
        values = zeros(n_points, n_runs)
        for (j, df) in enumerate(dfs)
            values[:, j] = df[1:n_points, col]
        end

        # Compute statistics
        result[!, "$(col)_mean"] = vec(mean(values, dims=2))
        result[!, "$(col)_std"] = vec(std(values, dims=2))
        result[!, "$(col)_sem"] = vec(std(values, dims=2) ./ sqrt(n_runs))
    end

    result[!, "n_runs"] = fill(n_runs, n_points)

    return result
end

"""
Aggregate with interpolation to common time points
"""
function aggregate_with_interpolation(dfs::Vector{DataFrame}, time_col::String,
                                      numeric_cols::Vector{String}, n_runs::Int)
    # Find common time range
    t_min = maximum(minimum(df[!, time_col]) for df in dfs)
    t_max = minimum(maximum(df[!, time_col]) for df in dfs)

    # Use time points from the file with most points in common range
    best_times = Float64[]
    for df in dfs
        times = df[!, time_col]
        valid_times = times[(times .>= t_min) .& (times .<= t_max)]
        if length(valid_times) > length(best_times)
            best_times = valid_times
        end
    end

    n_points = length(best_times)

    # Create result dataframe
    result = DataFrame()
    result[!, time_col] = best_times

    for col in numeric_cols
        # Interpolate values from all runs to common time points
        values = zeros(n_points, n_runs)
        for (j, df) in enumerate(dfs)
            t_orig = df[!, time_col]
            v_orig = df[!, col]

            # Linear interpolation
            for (i, t) in enumerate(best_times)
                idx = searchsortedlast(t_orig, t)
                if idx == 0
                    values[i, j] = v_orig[1]
                elseif idx == length(t_orig)
                    values[i, j] = v_orig[end]
                elseif t_orig[idx] == t
                    values[i, j] = v_orig[idx]
                else
                    # Linear interpolation
                    t1, t2 = t_orig[idx], t_orig[idx+1]
                    v1, v2 = v_orig[idx], v_orig[idx+1]
                    values[i, j] = v1 + (v2 - v1) * (t - t1) / (t2 - t1)
                end
            end
        end

        # Compute statistics
        result[!, "$(col)_mean"] = vec(mean(values, dims=2))
        result[!, "$(col)_std"] = vec(std(values, dims=2))
        result[!, "$(col)_sem"] = vec(std(values, dims=2) ./ sqrt(n_runs))
    end

    result[!, "n_runs"] = fill(n_runs, n_points)

    return result
end

"""
Save aggregated results to CSV
"""
function save_aggregated(df::DataFrame, output_path::String)
    # Create output directory if needed
    output_dir = dirname(output_path)
    if !isempty(output_dir) && !isdir(output_dir)
        mkpath(output_dir)
    end

    CSV.write(output_path, df)
    println("Saved aggregated data to: $output_path")
end

"""
Main function
"""
function main()
    args = parse_commandline()

    # Validate arguments
    if isempty(args["pattern"]) && isempty(args["base"])
        error("Must specify either --pattern or --base")
    end

    if !isempty(args["pattern"]) && !isempty(args["base"])
        error("Cannot specify both --pattern and --base")
    end

    interpolate = args["interpolate"]

    if !isempty(args["pattern"])
        # Single pattern mode
        pattern = args["pattern"]
        files = find_matching_files(pattern)

        println("Found $(length(files)) files matching pattern:")
        for f in files
            println("  - $f")
        end

        # Determine output path
        output = args["output"]
        if isempty(output)
            # Generate output name from pattern
            base_pattern = replace(basename(pattern), "*" => "aggregated")
            output = joinpath(dirname(pattern), base_pattern)
        end

        # Read and aggregate
        dfs = read_csv_files(files)
        result = aggregate_dataframes(dfs; interpolate=interpolate)
        save_aggregated(result, output)

    else
        # Base name mode - aggregate all metrics
        base = args["base"]
        phase = args["phase"]
        input_dir = args["input-dir"]
        output_dir = args["output-dir"]

        metrics = find_metric_files(base, phase, input_dir)

        println("Found metrics for base='$base', phase='$phase':")
        for (metric, files) in metrics
            println("  $metric: $(length(files)) files")
        end
        println()

        # Create output directory
        if !isdir(output_dir)
            mkpath(output_dir)
        end

        # Aggregate each metric type
        for (metric, files) in metrics
            println("Aggregating $metric from $(length(files)) runs...")
            for f in files
                println("  - $(basename(f))")
            end

            dfs = read_csv_files(files)
            result = aggregate_dataframes(dfs; interpolate=interpolate)

            output_file = joinpath(output_dir, "$(base)_$(phase)_$(metric)_aggregated.csv")
            save_aggregated(result, output_file)
            println()
        end
    end

    println("Aggregation complete!")
end

# Run main
main()
