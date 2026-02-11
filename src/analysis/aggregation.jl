# =============================================================================
# aggregation.jl - Reusable functions for aggregating sweep data
# =============================================================================
#
# This module provides functions for:
# - Parsing sweep filenames to extract parameters
# - Time-series aggregation across replicas
# - Steady-state summary statistics
# - Grouping files by parameter values
#
# =============================================================================

using CSV
using DataFrames
using Statistics

export parse_sweep_filename, extract_run_index
export detect_time_column, get_value_columns
export aggregate_timeseries, aggregate_summary
export group_files_by_params, compute_tail_statistics

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

# Example
```julia
params = parse_sweep_filename("single_N200_Nact50_kang1.0_fact5.0_run0_active_rg.csv")
# (system=:single, n_monomers=200, n_active=50, kangle=1.0, fact=5.0, simid="run0", phase=:active, metric="rg")
```
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

# =============================================================================
# Column Detection
# =============================================================================

"""
    detect_time_column(df::DataFrame)

Detect the time/index column in a DataFrame.
Checks for common time column names: time, lag_time, t, step, s, sample.
Returns the first matching column name, or the first column if none match.
"""
function detect_time_column(df::DataFrame)
    possible_names = ["time", "lag_time", "t", "step", "s", "sample"]
    for name in possible_names
        if name in names(df)
            return name
        end
    end
    return names(df)[1]
end

"""
    get_value_columns(df::DataFrame, time_col::String)

Get all numeric value columns (excluding the time column).
"""
function get_value_columns(df::DataFrame, time_col::String)
    return [name for name in names(df)
            if name != time_col && eltype(df[!, name]) <: Union{Number, Missing}]
end

# =============================================================================
# File Grouping
# =============================================================================

"""
    group_files_by_params(files::Vector{String}, group_by::Vector{Symbol})

Group files by specified parameter values.

# Arguments
- `files`: Vector of file paths
- `group_by`: Vector of parameter symbols to group by (e.g., [:n_active, :kangle])

# Returns
Dict{NamedTuple => Vector{String}} mapping parameter combinations to file lists.

# Example
```julia
files = ["single_N200_Nact0_kang1.0_fact5.0_run0_active_rg.csv",
         "single_N200_Nact0_kang1.0_fact5.0_run1_active_rg.csv",
         "single_N200_Nact50_kang1.0_fact5.0_run0_active_rg.csv"]
groups = group_files_by_params(files, [:n_active])
# Dict((n_active=0,) => ["..._Nact0_..._run0_...", "..._Nact0_..._run1_..."],
#      (n_active=50,) => ["..._Nact50_..._run0_..."])
```
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

        # Create grouping key
        values = [getfield(params, sym) for sym in group_by]
        key = NamedTuple{Tuple(group_by)}(values)

        if !haskey(groups, key)
            groups[key] = String[]
        end
        push!(groups[key], file)
    end

    return groups
end

# =============================================================================
# Time-Series Aggregation
# =============================================================================

"""
    aggregate_timeseries(files::Vector{String}; interpolate::Bool=false)
    aggregate_timeseries(dfs::Vector{DataFrame}; interpolate::Bool=false)

Aggregate multiple time-series into a single DataFrame with statistics.

For each value column in the input, computes:
- `{col}_mean`: Mean across all files at each time point
- `{col}_std`: Standard deviation
- `{col}_sem`: Standard error of the mean

Also adds `n_runs` column with the number of files aggregated.

# Arguments
- `files` or `dfs`: Input file paths or DataFrames
- `interpolate`: If true, interpolate to common time points when times don't match exactly

# Example
```julia
result = aggregate_timeseries(["run0_rg.csv", "run1_rg.csv", "run2_rg.csv"])
# Returns DataFrame with columns: time, Rg_mean, Rg_std, Rg_sem, n_runs
```
"""
function aggregate_timeseries(files::Vector{String}; interpolate::Bool=false)
    if isempty(files)
        error("No files to aggregate")
    end
    dfs = [CSV.read(f, DataFrame) for f in files]
    return aggregate_timeseries(dfs; interpolate=interpolate)
end

function aggregate_timeseries(dfs::Vector{DataFrame}; interpolate::Bool=false)
    if isempty(dfs)
        error("No DataFrames to aggregate")
    end

    n_runs = length(dfs)
    time_col = detect_time_column(first(dfs))
    value_cols = get_value_columns(first(dfs), time_col)

    if interpolate
        return _aggregate_with_interpolation(dfs, time_col, value_cols, n_runs)
    else
        return _aggregate_exact_match(dfs, time_col, value_cols, n_runs)
    end
end

function _aggregate_exact_match(dfs::Vector{DataFrame}, time_col::String,
                                 value_cols::Vector{String}, n_runs::Int)
    n_points = minimum(nrow.(dfs))

    max_rows = maximum(nrow.(dfs))
    if n_points != max_rows
        @warn "Files have different lengths ($n_points to $max_rows rows), using minimum"
    end

    result = DataFrame()
    result[!, time_col] = dfs[1][1:n_points, time_col]

    for col in value_cols
        values = zeros(n_points, n_runs)
        for (j, df) in enumerate(dfs)
            col_data = df[1:n_points, col]
            values[:, j] = coalesce.(col_data, NaN)
        end

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

function _aggregate_with_interpolation(dfs::Vector{DataFrame}, time_col::String,
                                        value_cols::Vector{String}, n_runs::Int)
    t_min = maximum(minimum(df[!, time_col]) for df in dfs)
    t_max = minimum(maximum(df[!, time_col]) for df in dfs)

    best_times = Float64[]
    for df in dfs
        times = df[!, time_col]
        valid_times = times[(times .>= t_min) .& (times .<= t_max)]
        if length(valid_times) > length(best_times)
            best_times = collect(valid_times)
        end
    end

    n_points = length(best_times)
    result = DataFrame()
    result[!, time_col] = best_times

    for col in value_cols
        values = zeros(n_points, n_runs)
        for (j, df) in enumerate(dfs)
            t_orig = df[!, time_col]
            v_orig = coalesce.(df[!, col], NaN)

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
    compute_tail_statistics(df::DataFrame; tail_percent::Float64=20.0)
    compute_tail_statistics(df::DataFrame, time_col::String, value_cols::Vector{String};
                            tail_percent::Float64=20.0)

Compute statistics from the tail portion of a time series.

# Arguments
- `df`: Input DataFrame
- `tail_percent`: Percentage of data to use from the end (default: 20%)

# Returns
Dict{Symbol, Float64} with `{col}_mean` and `{col}_std` for each value column.
"""
function compute_tail_statistics(df::DataFrame; tail_percent::Float64=20.0)
    time_col = detect_time_column(df)
    value_cols = get_value_columns(df, time_col)
    return compute_tail_statistics(df, time_col, value_cols; tail_percent=tail_percent)
end

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

Generate a summary table with steady-state statistics.

For each file, extracts tail statistics and file parameters. If `group_by` is specified,
aggregates statistics across replicas within each parameter group.

# Arguments
- `files`: Vector of file paths
- `tail_percent`: Percentage of data for tail statistics (default: 20%)
- `group_by`: Parameters to group by when aggregating across replicas

# Returns
DataFrame with one row per file (or per parameter group if grouping).

# Example
```julia
# Get per-file statistics
summary = aggregate_summary(files; tail_percent=20.0)

# Aggregate across replicas, one row per n_active value
summary = aggregate_summary(files; tail_percent=20.0, group_by=[:n_active])
```
"""
function aggregate_summary(files::Vector{String}; tail_percent::Float64=20.0,
                            group_by::Vector{Symbol}=Symbol[])
    if isempty(files)
        error("No files to aggregate")
    end

    file_stats = Dict{Symbol, Any}[]

    for file in files
        df = CSV.read(file, DataFrame)
        params = parse_sweep_filename(file)

        if params === nothing
            @warn "Could not parse filename: $(basename(file)), skipping"
            continue
        end

        tail_stats = compute_tail_statistics(df; tail_percent=tail_percent)

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

    summary_df = DataFrame(file_stats)

    if !isempty(group_by)
        summary_df = _aggregate_across_replicas(summary_df, group_by)
    end

    return summary_df
end

function _aggregate_across_replicas(df::DataFrame, group_by::Vector{Symbol})
    mean_cols = [name for name in names(df) if endswith(string(name), "_mean")]
    metric_bases = [replace(string(col), "_mean" => "") for col in mean_cols]

    grouped = groupby(df, group_by)

    results = DataFrame[]
    for g in grouped
        row = Dict{Symbol, Any}()

        for sym in group_by
            row[sym] = first(g[!, sym])
        end

        row[:n_runs] = nrow(g)

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
