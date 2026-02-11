#!/usr/bin/env julia
# =============================================================================
# aggregate_sweep.jl - Aggregate simulation data from parameter sweeps
# =============================================================================
#
# A unified tool for aggregating simulation data with three modes:
#
# Usage:
#   # Mode 1: Time-series aggregation (default)
#   julia --project=. scripts/aggregate_sweep.jl \
#       --input-dir _data/csv \
#       --pattern "single_N200_Nact*_run*_active_rg.csv" \
#       --output-dir _data/aggregated \
#       --mode timeseries
#
#   # Mode 2: Steady-state summary
#   julia --project=. scripts/aggregate_sweep.jl \
#       --input-dir _data/csv \
#       --pattern "single_N200_Nact*_run*_active_rg.csv" \
#       --output _data/aggregated/rg_summary.csv \
#       --mode summary \
#       --tail-percent 20
#
#   # Mode 3: Per-parameter grouping
#   julia --project=. scripts/aggregate_sweep.jl \
#       --input-dir _data/csv \
#       --pattern "single_N200_Nact*_run*_active_rg.csv" \
#       --output-dir _data/aggregated \
#       --mode timeseries \
#       --group-by n_active
#
# =============================================================================

using ArgParse
using CSV
using DataFrames
using Statistics
using Glob

# =============================================================================
# Filename Parsing
# =============================================================================

"""
    parse_sweep_filename(filename::String)

Extract parameters from a sweep filename.

Supports patterns like:
- `single_N{n}_Nact{nact}_kang{k}_fact{f}_{simid}_{phase}_{metric}.csv`
- `double_N1{n1}_N2{n2}_Nact1{na1}_Nact2{na2}_kang{k}_fact{f}_{simid}_{phase}_{metric}.csv`

Returns a NamedTuple with extracted parameters, or nothing if pattern doesn't match.
"""
function parse_sweep_filename(filename::String)
    basename_str = basename(filename)

    # Try single ring pattern first
    single_pattern = r"^single_N(\d+)_Nact(\d+)_kang([\d.]+)_fact([\d.]+)_([^_]+)_(\w+)_(\w+)\.csv$"
    m = match(single_pattern, basename_str)
    if m !== nothing
        return (
            system = :single,
            n_monomers = parse(Int, m[1]),
            n_active = parse(Int, m[2]),
            kangle = parse(Float64, m[3]),
            fact = parse(Float64, m[4]),
            simid = String(m[5]),
            phase = Symbol(m[6]),
            metric = String(m[7])
        )
    end

    # Try double ring pattern
    double_pattern = r"^double_N1(\d+)_N2(\d+)_Nact1(\d+)_Nact2(\d+)_kang([\d.]+)_fact([\d.]+)_([^_]+)_(\w+)_(\w+)\.csv$"
    m = match(double_pattern, basename_str)
    if m !== nothing
        return (
            system = :double,
            n_monomers_1 = parse(Int, m[1]),
            n_monomers_2 = parse(Int, m[2]),
            n_active_1 = parse(Int, m[3]),
            n_active_2 = parse(Int, m[4]),
            kangle = parse(Float64, m[5]),
            fact = parse(Float64, m[6]),
            simid = String(m[7]),
            phase = Symbol(m[8]),
            metric = String(m[9])
        )
    end

    return nothing
end

"""
    extract_run_index(simid::String)

Extract run index from simid (e.g., "run0" -> 0, "run10" -> 10).
Returns -1 if not a run-indexed simid.
"""
function extract_run_index(simid::String)
    m = match(r"^run(\d+)$", simid)
    return m !== nothing ? parse(Int, m[1]) : -1
end

"""
    get_grouping_key(params::NamedTuple, group_by::Vector{Symbol})

Create a grouping key from parameters based on specified fields.
"""
function get_grouping_key(params::NamedTuple, group_by::Vector{Symbol})
    values = [getfield(params, sym) for sym in group_by]
    return NamedTuple{Tuple(group_by)}(values)
end

"""
    group_files_by_params(files::Vector{String}, group_by::Vector{Symbol})

Group files by specified parameter values.
Returns Dict{NamedTuple => Vector{String}}.
"""
function group_files_by_params(files::Vector{String}, group_by::Vector{Symbol})
    groups = Dict{Any, Vector{String}}()

    for file in files
        params = parse_sweep_filename(file)
        if params === nothing
            @warn "Could not parse filename: $(basename(file)), skipping"
            continue
        end

        # Check that all group_by fields exist
        valid = true
        for sym in group_by
            if !haskey(params, sym)
                @warn "Parameter $sym not found in $(basename(file)), skipping"
                valid = false
                break
            end
        end
        if !valid
            continue
        end

        key = get_grouping_key(params, group_by)
        if !haskey(groups, key)
            groups[key] = String[]
        end
        push!(groups[key], file)
    end

    return groups
end

# =============================================================================
# Time Column Detection
# =============================================================================

"""
    detect_time_column(df::DataFrame)

Detect the time/index column in a DataFrame.
Returns the column name.
"""
function detect_time_column(df::DataFrame)
    possible_names = ["time", "lag_time", "t", "step", "s", "sample"]
    for name in possible_names
        if name in names(df)
            return name
        end
    end
    # Default to first column
    return names(df)[1]
end

"""
    get_value_columns(df::DataFrame, time_col::String)

Get all numeric value columns (excluding time column).
"""
function get_value_columns(df::DataFrame, time_col::String)
    return [name for name in names(df)
            if name != time_col && eltype(df[!, name]) <: Union{Number, Missing}]
end

# =============================================================================
# Time-Series Aggregation
# =============================================================================

"""
    aggregate_timeseries(files::Vector{String}; interpolate::Bool=false)

Aggregate multiple CSV files into a single time-series with statistics.

Returns a DataFrame with columns:
- time column (preserved from input)
- {col}_mean, {col}_std, {col}_sem for each value column
- n_runs: number of files aggregated
"""
function aggregate_timeseries(files::Vector{String}; interpolate::Bool=false)
    if isempty(files)
        error("No files to aggregate")
    end

    dfs = [CSV.read(f, DataFrame) for f in files]
    n_runs = length(dfs)

    # Detect time column and value columns from first file
    time_col = detect_time_column(first(dfs))
    value_cols = get_value_columns(first(dfs), time_col)

    if interpolate
        return aggregate_with_interpolation(dfs, time_col, value_cols, n_runs)
    else
        return aggregate_exact_match(dfs, time_col, value_cols, n_runs)
    end
end

"""
    aggregate_exact_match(dfs, time_col, value_cols, n_runs)

Aggregate assuming exact time point matching across files.
"""
function aggregate_exact_match(dfs::Vector{DataFrame}, time_col::String,
                                value_cols::Vector{String}, n_runs::Int)
    # Find minimum number of rows
    n_points = minimum(nrow.(dfs))

    # Warn if files have different lengths
    max_rows = maximum(nrow.(dfs))
    if n_points != max_rows
        @warn "Files have different lengths ($n_points to $max_rows rows), using minimum"
    end

    # Create result DataFrame
    result = DataFrame()
    result[!, time_col] = dfs[1][1:n_points, time_col]

    for col in value_cols
        # Collect values from all runs
        values = zeros(n_points, n_runs)
        for (j, df) in enumerate(dfs)
            col_data = df[1:n_points, col]
            # Handle missing values
            values[:, j] = coalesce.(col_data, NaN)
        end

        # Compute statistics (ignoring NaN)
        result[!, "$(col)_mean"] = vec(mapslices(x -> mean(filter(!isnan, x)), values, dims=2))
        result[!, "$(col)_std"] = vec(mapslices(x -> std(filter(!isnan, x)), values, dims=2))
        result[!, "$(col)_sem"] = vec(mapslices(x -> begin
            valid = filter(!isnan, x)
            isempty(valid) ? NaN : std(valid) / sqrt(length(valid))
        end, values, dims=2))
    end

    result[!, :n_runs] .= n_runs

    return result
end

"""
    aggregate_with_interpolation(dfs, time_col, value_cols, n_runs)

Aggregate with interpolation to common time points.
"""
function aggregate_with_interpolation(dfs::Vector{DataFrame}, time_col::String,
                                       value_cols::Vector{String}, n_runs::Int)
    # Find common time range
    t_min = maximum(minimum(df[!, time_col]) for df in dfs)
    t_max = minimum(maximum(df[!, time_col]) for df in dfs)

    # Use time points from the file with most points in common range
    best_times = Float64[]
    for df in dfs
        times = df[!, time_col]
        valid_times = times[(times .>= t_min) .& (times .<= t_max)]
        if length(valid_times) > length(best_times)
            best_times = collect(valid_times)
        end
    end

    n_points = length(best_times)

    # Create result DataFrame
    result = DataFrame()
    result[!, time_col] = best_times

    for col in value_cols
        values = zeros(n_points, n_runs)
        for (j, df) in enumerate(dfs)
            t_orig = df[!, time_col]
            v_orig = coalesce.(df[!, col], NaN)

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
                    t1, t2 = t_orig[idx], t_orig[idx+1]
                    v1, v2 = v_orig[idx], v_orig[idx+1]
                    values[i, j] = v1 + (v2 - v1) * (t - t1) / (t2 - t1)
                end
            end
        end

        # Compute statistics
        result[!, "$(col)_mean"] = vec(mapslices(x -> mean(filter(!isnan, x)), values, dims=2))
        result[!, "$(col)_std"] = vec(mapslices(x -> std(filter(!isnan, x)), values, dims=2))
        result[!, "$(col)_sem"] = vec(mapslices(x -> begin
            valid = filter(!isnan, x)
            isempty(valid) ? NaN : std(valid) / sqrt(length(valid))
        end, values, dims=2))
    end

    result[!, :n_runs] .= n_runs

    return result
end

# =============================================================================
# Steady-State Summary
# =============================================================================

"""
    compute_tail_statistics(df::DataFrame, time_col::String, value_cols::Vector{String};
                            tail_percent::Float64=20.0)

Compute statistics from the tail portion of a time series.
"""
function compute_tail_statistics(df::DataFrame, time_col::String, value_cols::Vector{String};
                                  tail_percent::Float64=20.0)
    n = nrow(df)
    tail_start = max(1, n - round(Int, n * tail_percent / 100) + 1)
    tail_data = df[tail_start:end, :]

    stats = Dict{Symbol, Float64}()
    for col in value_cols
        col_data = coalesce.(tail_data[!, col], NaN)
        valid = filter(!isnan, col_data)
        if !isempty(valid)
            stats[Symbol("$(col)_mean")] = mean(valid)
            stats[Symbol("$(col)_std")] = std(valid)
        else
            stats[Symbol("$(col)_mean")] = NaN
            stats[Symbol("$(col)_std")] = NaN
        end
    end

    return stats
end

"""
    aggregate_summary(files::Vector{String}; tail_percent::Float64=20.0,
                      group_by::Vector{Symbol}=Symbol[])

Generate a summary table with steady-state statistics for each file,
optionally aggregated across replicas within parameter groups.

Returns a DataFrame with one row per parameter combination.
"""
function aggregate_summary(files::Vector{String}; tail_percent::Float64=20.0,
                            group_by::Vector{Symbol}=Symbol[])
    if isempty(files)
        error("No files to aggregate")
    end

    # First pass: collect per-file statistics
    file_stats = Dict{Symbol, Any}[]

    for file in files
        df = CSV.read(file, DataFrame)
        params = parse_sweep_filename(file)

        if params === nothing
            @warn "Could not parse filename: $(basename(file)), skipping"
            continue
        end

        time_col = detect_time_column(df)
        value_cols = get_value_columns(df, time_col)

        tail_stats = compute_tail_statistics(df, time_col, value_cols; tail_percent=tail_percent)

        # Combine params and stats
        row = Dict{Symbol, Any}()
        for (k, v) in pairs(params)
            row[k] = v
        end
        for (k, v) in tail_stats
            row[k] = v
        end
        push!(file_stats, row)
    end

    if isempty(file_stats)
        error("No valid files found")
    end

    # Convert to DataFrame
    summary_df = DataFrame(file_stats)

    # If group_by is specified, aggregate across replicas
    if !isempty(group_by)
        summary_df = aggregate_across_replicas(summary_df, group_by)
    end

    return summary_df
end

"""
    aggregate_across_replicas(df::DataFrame, group_by::Vector{Symbol})

Aggregate statistics across replicas within each parameter group.
"""
function aggregate_across_replicas(df::DataFrame, group_by::Vector{Symbol})
    # Find all _mean and _std columns (these are per-file tail statistics)
    mean_cols = [name for name in names(df) if endswith(string(name), "_mean")]

    # Determine the base metric names
    metric_bases = [replace(string(col), "_mean" => "") for col in mean_cols]

    # Group by parameters
    grouped = groupby(df, group_by)

    results = DataFrame[]
    for g in grouped
        row = Dict{Symbol, Any}()

        # Copy grouping keys
        for sym in group_by
            row[sym] = first(g[!, sym])
        end

        # Count replicas
        row[:n_runs] = nrow(g)

        # Aggregate each metric
        for base in metric_bases
            mean_col = Symbol("$(base)_mean")
            if mean_col in propertynames(g)
                values = g[!, mean_col]
                valid = filter(!isnan, values)
                if !isempty(valid)
                    row[mean_col] = mean(valid)
                    row[Symbol("$(base)_std")] = std(valid)
                    row[Symbol("$(base)_sem")] = std(valid) / sqrt(length(valid))
                end
            end
        end

        push!(results, DataFrame([row]))
    end

    return vcat(results...)
end

# =============================================================================
# Output Naming
# =============================================================================

"""
    generate_output_name(pattern::String, mode::Symbol; group_key=nothing, metric::String="")

Generate output filename based on input pattern and mode.
"""
function generate_output_name(pattern::String, mode::Symbol; group_key=nothing, metric::String="")
    # Clean up pattern to create base name
    base = basename(pattern)
    base = replace(base, r"\*" => "")
    base = replace(base, r"_run\d*" => "")
    base = replace(base, r"__+" => "_")  # Remove double underscores
    base = replace(base, r"^_|_$" => "") # Remove leading/trailing underscores
    base = replace(base, ".csv" => "")

    suffix = if mode == :summary
        "_summary"
    else
        "_aggregated"
    end

    if group_key !== nothing
        group_str = join(["$(k)$(v)" for (k, v) in pairs(group_key)], "_")
        return "$(base)_$(group_str)$(suffix).csv"
    else
        return "$(base)$(suffix).csv"
    end
end

# =============================================================================
# File Discovery
# =============================================================================

"""
    find_matching_files(input_dir::String, pattern::String; metric::String="")

Find all files matching the pattern, optionally filtered by metric.
If the pattern already contains a directory path, input_dir is ignored.
"""
function find_matching_files(input_dir::String, pattern::String; metric::String="")
    # Check if pattern already contains a directory path
    full_pattern = if occursin("/", pattern) || occursin("\\", pattern)
        pattern  # Pattern already has a path, use as-is
    else
        joinpath(input_dir, pattern)
    end
    files = glob(full_pattern)

    if isempty(files)
        error("No files found matching pattern: $full_pattern")
    end

    # Filter by metric if specified
    if !isempty(metric)
        files = filter(f -> begin
            params = parse_sweep_filename(f)
            params !== nothing && params.metric == metric
        end, files)

        if isempty(files)
            error("No files found with metric '$metric'")
        end
    end

    return sort(files)
end

# =============================================================================
# CLI Interface
# =============================================================================

function parse_commandline()
    s = ArgParseSettings(
        description = "Aggregate simulation data from parameter sweeps with multiple replicas",
        prog = "aggregate_sweep.jl"
    )

    @add_arg_table! s begin
        "--input-dir"
            help = "Directory containing CSV files"
            arg_type = String
            default = "_data/csv"
        "--pattern"
            help = "Glob pattern to match files (required)"
            arg_type = String
            required = true
        "--output-dir"
            help = "Output directory for aggregated files"
            arg_type = String
            default = "_data/aggregated"
        "--output"
            help = "Single output file path (for summary mode without grouping)"
            arg_type = String
            default = ""
        "--mode"
            help = "Aggregation mode: timeseries, summary, or both"
            arg_type = String
            default = "timeseries"
        "--group-by"
            help = "Parameter(s) to group by (comma-separated, e.g., 'n_active' or 'n_active,kangle')"
            arg_type = String
            default = ""
        "--tail-percent"
            help = "Percentage of data for tail statistics (summary mode)"
            arg_type = Float64
            default = 20.0
        "--interpolate"
            help = "Interpolate to common time points"
            action = :store_true
        "--metric"
            help = "Filter to specific metric (e.g., 'rg', 'msd')"
            arg_type = String
            default = ""
        "--verbose"
            help = "Print detailed progress information"
            action = :store_true
    end

    return parse_args(s)
end

function save_result(df::DataFrame, output_path::String)
    output_dir = dirname(output_path)
    if !isempty(output_dir) && !isdir(output_dir)
        mkpath(output_dir)
    end
    CSV.write(output_path, df)
    println("  Saved: $output_path")
end

function main()
    args = parse_commandline()

    # Parse arguments
    input_dir = args["input-dir"]
    pattern = args["pattern"]
    output_dir = args["output-dir"]
    output = args["output"]
    mode = Symbol(args["mode"])
    tail_percent = args["tail-percent"]
    interpolate = args["interpolate"]
    metric = args["metric"]
    verbose = args["verbose"]

    # Parse group-by parameters
    group_by_str = args["group-by"]
    group_by = isempty(group_by_str) ? Symbol[] : Symbol.(split(group_by_str, ","))

    # Validate mode
    if mode ∉ [:timeseries, :summary, :both]
        error("Invalid mode: $mode. Must be 'timeseries', 'summary', or 'both'")
    end

    # Find matching files
    # If pattern contains a path, use it directly; otherwise prepend input_dir
    search_pattern = if occursin("/", pattern) || occursin("\\", pattern)
        pattern
    else
        joinpath(input_dir, pattern)
    end
    println("Searching for files matching: $search_pattern")
    files = find_matching_files(input_dir, pattern; metric=metric)
    println("Found $(length(files)) files")

    if verbose
        for f in files
            println("  - $(basename(f))")
        end
    end

    # Create output directory
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Process based on mode and grouping
    if isempty(group_by)
        # No grouping - aggregate all files together
        println("\nAggregating all $(length(files)) files...")

        if mode == :timeseries || mode == :both
            println("Computing time-series aggregation...")
            result = aggregate_timeseries(files; interpolate=interpolate)
            out_path = !isempty(output) && mode == :timeseries ?
                       output : joinpath(output_dir, generate_output_name(pattern, :timeseries))
            save_result(result, out_path)
        end

        if mode == :summary || mode == :both
            println("Computing summary statistics (tail $(tail_percent)%)...")
            result = aggregate_summary(files; tail_percent=tail_percent)
            out_path = !isempty(output) && mode == :summary ?
                       output : joinpath(output_dir, generate_output_name(pattern, :summary))
            save_result(result, out_path)
        end
    else
        # Group files by parameter values
        println("\nGrouping files by: $(join(string.(group_by), ", "))")
        groups = group_files_by_params(files, group_by)
        println("Found $(length(groups)) groups")

        # For summary mode with grouping, collect all results to combine into single file
        summary_results = DataFrame[]

        for (key, group_files) in sort(collect(groups), by=x -> x[1])
            key_str = join(["$(k)=$(v)" for (k, v) in pairs(key)], ", ")
            println("\n  Group [$key_str]: $(length(group_files)) files")

            if verbose
                for f in group_files
                    println("    - $(basename(f))")
                end
            end

            if mode == :timeseries || mode == :both
                result = aggregate_timeseries(group_files; interpolate=interpolate)
                out_path = joinpath(output_dir, generate_output_name(pattern, :timeseries; group_key=key))
                save_result(result, out_path)
            end

            if mode == :summary || mode == :both
                result = aggregate_summary(group_files; tail_percent=tail_percent, group_by=group_by)
                push!(summary_results, result)
            end
        end

        # Save combined summary if in summary mode
        if (mode == :summary || mode == :both) && !isempty(summary_results)
            combined_summary = vcat(summary_results...)
            out_path = !isempty(output) ? output : joinpath(output_dir, generate_output_name(pattern, :summary))
            save_result(combined_summary, out_path)
        end
    end

    println("\nAggregation complete!")
end

# Run main
main()
